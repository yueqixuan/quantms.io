"""LazyQuery builder — immutable, chainable SQL query construction."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from qpx.core.engine import DuckDBEngine


def _escape_sql_string(value: str) -> str:
    """Escape single quotes in a SQL string literal to prevent injection."""
    return value.replace("'", "''")


class LazyQuery:
    """
    Immutable lazy SQL query builder.

    Every filter/select/join returns a NEW LazyQuery — the original is
    unchanged. Nothing hits DuckDB until materialize() is called.

    This is the query engine shared by ALL data structures. Optimizations
    (predicate pushdown hints, query caching, etc.) go here.
    """

    def __init__(self, engine: "DuckDBEngine", table_name: str, sql: str | None = None):
        self._engine = engine
        self._table_name = table_name
        self._sql = sql  # None means "SELECT * FROM table_name"

    @property
    def source(self) -> str:
        """SQL fragment for use in FROM clause."""
        if self._sql:
            return f"({self._sql})"
        return self._table_name

    def select(self, *columns: str) -> LazyQuery:
        cols = ", ".join(columns)
        return LazyQuery(self._engine, self._table_name, f"SELECT {cols} FROM {self.source}")

    def filter(self, condition: str) -> LazyQuery:
        return LazyQuery(
            self._engine,
            self._table_name,
            f"SELECT * FROM {self.source} WHERE {condition}",
        )

    def limit(self, n: int) -> LazyQuery:
        return LazyQuery(self._engine, self._table_name, f"SELECT * FROM {self.source} LIMIT {n}")

    def join(self, other: LazyQuery, on: str) -> LazyQuery:
        sql = f"SELECT * FROM {self.source} a JOIN {other.source} b ON a.{on} = b.{on}"
        return LazyQuery(self._engine, self._table_name, sql)

    def order_by(self, *columns: str, desc: bool = False) -> LazyQuery:
        direction = "DESC" if desc else "ASC"
        cols = ", ".join(f"{c} {direction}" for c in columns)
        return LazyQuery(
            self._engine,
            self._table_name,
            f"SELECT * FROM {self.source} ORDER BY {cols}",
        )

    def group_by(self, *columns: str, agg: str = "COUNT(*) AS count") -> LazyQuery:
        cols = ", ".join(columns)
        return LazyQuery(
            self._engine,
            self._table_name,
            f"SELECT {cols}, {agg} FROM {self.source} GROUP BY {cols}",
        )

    def distinct_values(self, column: str) -> list:
        """Return distinct values for a column."""
        result = self._engine.execute(f"SELECT DISTINCT {column} FROM {self.source} ORDER BY {column}")
        return [row[0] for row in result.fetchall()]

    def build_sql(self) -> str:
        if self._sql:
            return self._sql
        return f"SELECT * FROM {self._table_name}"

    def execute(self):
        """Execute the query and return raw DuckDB result."""
        return self._engine.execute(self.build_sql())

    def count(self) -> int:
        row = self._engine.execute(f"SELECT COUNT(*) FROM {self.source}").fetchone()
        return row[0]
