"""Module for storing type definitions used in the FFCx code base."""

from typing import Literal

entity_types = Literal["cell", "facet", "ridge", "peak"]

supported_integral_types = Literal["cell", "interior_facet", "exterior_facet", "ridge", "peak"]
