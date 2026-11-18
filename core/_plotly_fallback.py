"""
Lightweight Plotly Figure/Trace fallback for environments where `plotly` is not installed.
When `plotly` is installed, MolecuLens imports and uses `plotly.graph_objects` directly;
this module is loaded ONLY in `except ImportError:` fallback paths so backend engines and
test suites can execute without optional charting packages.
"""

from __future__ import annotations

import sys
import types
from typing import Any


class _DictAttr(dict):
    """Dict that also supports attribute read/write access like Plotly layout objects."""

    def __getattr__(self, name: str) -> Any:
        try:
            val = self[name]
            if isinstance(val, dict) and not isinstance(val, _DictAttr):
                val = _DictAttr(val)
                self[name] = val
            return val
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __setattr__(self, name: str, value: Any) -> None:
        self[name] = _DictAttr(value) if isinstance(value, dict) and not isinstance(value, _DictAttr) else value


class _Trace:
    def __init__(self, trace_type: str = "scatter", **kwargs: Any):
        self.type = trace_type
        for k, v in kwargs.items():
            setattr(self, k, v)

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


class Scatter(_Trace):
    def __init__(self, **kwargs: Any):
        super().__init__("scatter", **kwargs)


class Scatter3d(_Trace):
    def __init__(self, **kwargs: Any):
        super().__init__("scatter3d", **kwargs)

class Bar(_Trace):
    def __init__(self, **kwargs: Any):
        super().__init__("bar", **kwargs)


class Figure:
    def __init__(self, data: Any = None, layout: Any = None):
        self.data: list[Any] = list(data) if data else []
        self.layout: _DictAttr = _DictAttr(
            layout
            or {
                "title": _DictAttr({"text": ""}),
                "xaxis": _DictAttr({"showticklabels": True}),
                "yaxis": _DictAttr({"showticklabels": True}),
                "annotations": [],
                "shapes": [],
            }
        )

    def add_trace(self, trace: Any) -> Figure:
        self.data.append(trace)
        return self

    def add_shape(self, **kwargs: Any) -> Figure:
        self.layout.setdefault("shapes", []).append(dict(kwargs))
        return self

    def add_annotation(self, **kwargs: Any) -> Figure:
        self.layout.setdefault("annotations", []).append(dict(kwargs))
        return self

    def add_hline(self, **kwargs: Any) -> Figure:
        self.layout.setdefault("shapes", []).append({"type": "hline", **kwargs})
        return self

    def add_vline(self, **kwargs: Any) -> Figure:
        self.layout.setdefault("shapes", []).append({"type": "vline", **kwargs})
        return self

    def add_hrect(self, **kwargs: Any) -> Figure:
        self.layout.setdefault("shapes", []).append({"type": "hrect", **kwargs})
        return self

    def update_layout(self, **kwargs: Any) -> Figure:
        for k, v in kwargs.items():
            if isinstance(v, dict):
                existing = self.layout.get(k)
                if isinstance(existing, dict):
                    existing.update(v)
                    self.layout[k] = _DictAttr(existing)
                else:
                    self.layout[k] = _DictAttr(v)
            else:
                self.layout[k] = v
        return self

    def to_dict(self) -> dict[str, Any]:
        return {
            "data": [t.to_dict() if hasattr(t, "to_dict") else dict(t.__dict__) for t in self.data],
            "layout": dict(self.layout),
        }


def register_plotly_fallback_in_sys_modules() -> None:
    if "plotly" in sys.modules:
        return
    plotly_mod = types.ModuleType("plotly")
    go_mod = types.ModuleType("plotly.graph_objects")
    go_mod.Figure = Figure
    go_mod.Scatter = Scatter
    go_mod.Bar = Bar
    go_mod.Scatter3d = Scatter3d
    plotly_mod.graph_objects = go_mod
    sys.modules["plotly"] = plotly_mod
    sys.modules["plotly.graph_objects"] = go_mod
