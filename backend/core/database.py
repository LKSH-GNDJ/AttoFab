"""
core/database.py - SQLAlchemy async ORM + auto-migrations.

"Auto-migration" here means `create_all()` on startup, deliberately - a
real migration tool (Alembic) is overkill for a single-table run-history
log at this project's current scope. If the schema grows enough to need
real migrations, that's the signal to introduce Alembic, not a reason to
add it preemptively.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
from typing import AsyncIterator

from sqlalchemy import DateTime, Integer, String, Text, select
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
DB_PATH = DATA_DIR / "attofab.db"

engine = create_async_engine(f"sqlite+aiosqlite:///{DB_PATH}", echo=False)
SessionLocal = async_sessionmaker(engine, expire_on_commit=False)


class Base(DeclarativeBase):
    pass


class Run(Base):
    """One simulation run: the recipe submitted, and a summary + full
    wafer JSON produced by the engine."""

    __tablename__ = "runs"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc)
    )
    nx: Mapped[int] = mapped_column(Integer)
    ny: Mapped[int] = mapped_column(Integer)
    recipe_json: Mapped[str] = mapped_column(Text)
    summary_json: Mapped[str] = mapped_column(Text)
    wafer_json: Mapped[str] = mapped_column(Text)

    def summary(self) -> dict:
        return {
            "id": self.id,
            "created_at": self.created_at.isoformat(),
            "nx": self.nx,
            "ny": self.ny,
            "summary": json.loads(self.summary_json),
        }


async def init_db() -> None:
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)


async def get_session() -> AsyncIterator[AsyncSession]:
    async with SessionLocal() as session:
        yield session


async def save_run(session: AsyncSession, *, nx: int, ny: int, recipe: dict, summary: dict, wafer: dict) -> Run:
    run = Run(
        nx=nx,
        ny=ny,
        recipe_json=json.dumps(recipe),
        summary_json=json.dumps(summary),
        wafer_json=json.dumps(wafer),
    )
    session.add(run)
    await session.commit()
    await session.refresh(run)
    return run


async def list_runs(session: AsyncSession, limit: int = 50) -> list[Run]:
    result = await session.execute(select(Run).order_by(Run.id.desc()).limit(limit))
    return list(result.scalars().all())


async def get_run(session: AsyncSession, run_id: int) -> Run | None:
    return await session.get(Run, run_id)
