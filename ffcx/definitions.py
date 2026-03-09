"""Module for storing type definitions used in the FFCx code base."""

from typing import Literal


entity_types = Literal["cell", "facet", "vertex", "ridge"]

supported_integral_types = Literal["cell", "interior_facet", "exterior_facet", "vertex", "ridge"]
