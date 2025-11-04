from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel


class ColumnSettings(BaseModel):
    use_as: str


class MapSettings(BaseModel):
    center: tuple[float, float]
    zoom_level: int


class Settings(BaseModel):
    master_columns: Optional[dict[str, ColumnSettings]] = None
    map: Optional[MapSettings] = None


def load_settings(settings_file: Path) -> Settings:
    """Load settings from a file."""
    data = yaml.safe_load(open(settings_file, "r"))
    return Settings(**data)


def get_master_columns_mapping(settings: Settings) -> dict[str, str]:
    if settings.master_columns is None:
        return {}
    return {col: column.use_as for col, column in settings.master_columns.items()}


def get_map_settings(
    settings: Settings,
) -> tuple[Optional[tuple[float, float]], Optional[int]]:
    if settings.map is None:
        return None, None
    return (settings.map.center, settings.map.zoom_level)
