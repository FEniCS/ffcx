"""Module for storing type definitions used in the FFCx code base."""

from functools import singledispatchmethod
from typing import Literal

import ufl
from ufl.corealg.dag_traverser import DAGTraverser

entity_types = Literal["cell", "facet", "vertex", "ridge"]

supported_integral_types = Literal["interior_facet", "exterior_facet", "ridge", "cell"]

class AveragedArgument(ufl.Argument):
    def __init__(self, u: ufl.Argument, average: Literal["cell", "facet"]):
        super().__init__(u.ufl_function_space(), u.number(), u.part())
        self._repr = "Averaged" + self._repr
        self._average = average

    def is_cellwise_constant(self)->bool:
        return True

    @property
    def average(self)->Literal["cell", "facet"]:
        return self._average


class IntermediateConstant(ufl.Constant):
    """A UFL constant used as a placeholder for arguments computed outside quadrature loops"""

    _original_shape: tuple[int, ...]

    def __init__(self, domain: ufl.domain.AbstractDomain, shape: tuple[int,...], count:int|None=None):
        super().__init__(domain, (), count)
        self._original_shape = shape

    @property
    def original_shape(self)-> tuple[int,...]:
        return self._original_shape


class AverageReplacer(DAGTraverser):
    """DAGTraverser to replaced averaged arguments with an argument in an
    intermediate space.
    """

    def __init__(
        self,
        compress: bool | None = True,
        visited_cache: dict[tuple, ufl.core.expr.Expr] | None = None,
        result_cache: dict[ufl.core.expr.Expr, ufl.core.expr.Expr] | None = None,
    ) -> None:
        """Initialise.

        Args:
            compress: If True, ``result_cache`` will be used.
            visited_cache: cache of intermediate results;
                expr -> r = self.process(expr, ...).
            result_cache: cache of result objects for memory reuse, r -> r.

        """
        self._sub_graphs: list[ufl.core.expr.Expr] = []
        super().__init__(compress=compress, visited_cache=visited_cache, result_cache=result_cache)

    @singledispatchmethod
    def process(
        self,
        o: ufl.core.expr.Expr,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Replace averaged arguments with intermediate space.

        Args:
            o: `ufl.core.expr.Expr` to be processed.
            reference_value: Whether `ReferenceValue` has been applied or not.
            reference_grad: Number of `ReferenceGrad`s that have been applied.
            restricted: '+', '-', or None.
        """
        return super().process(o)

    @process.register(ufl.averaging.CellAvg)
    def _(
        self,
        o: ufl.averaging.CellAvg,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Handle AverageOperator."""
        ops = o.ufl_operands
        assert len(ops) == 1, "Expected single operator in averaging"
        arguments = ufl.algorithms.extract_arguments(ops[0])
        replace_table = {}
        for argument in arguments:
            replace_table[argument] = AveragedArgument(argument, "cell")
        expr_domain = ufl.domain.extract_unique_domain(ops[0])
        cell_vol = ufl.CellVolume(expr_domain)
        cell_vol_lowered = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(
            cell_vol, (ufl.classes.Jacobian,)
        )
        self._sub_graphs.append(ops[0])
        assert expr_domain is not None
        constant = IntermediateConstant(expr_domain, shape=ops[0].ufl_shape)
        replaced_expr =  ufl.replace(ops[0], replace_table)
        new_expr = constant * replaced_expr / cell_vol_lowered
        assert isinstance(new_expr, ufl.core.expr.Expr)
        return new_expr

    @process.register(ufl.averaging.FacetAvg)
    def _(
        self,
        o: ufl.averaging.FacetAvg,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Handle AverageOperator."""
        ops = o.ufl_operands
        assert len(ops) == 1, "Expected single operator in averaging"
        arguments = ufl.algorithms.extract_arguments(ops[0])
        expr_domain = ufl.domain.extract_unique_domain(ops[0])
        assert expr_domain is not None
        replace_table = {}
        for argument in arguments:
            replace_table[argument] = AveragedArgument(argument, "facet")
        facet_area = ufl.FacetArea(expr_domain)
        facet_area_lowered = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(
            facet_area, (ufl.classes.Jacobian,)
        )
        self._sub_graphs.append(ops[0])
        constant = IntermediateConstant(expr_domain, shape=ops[0].ufl_shape)
        replaced_expr =  ufl.replace(ops[0], replace_table)
        new_expr = constant * replaced_expr/ facet_area_lowered
        assert isinstance(new_expr, ufl.core.expr.Expr)
        return new_expr

    @process.register(ufl.core.expr.Expr)
    def _(
        self,
        o: ufl.Argument,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Handle anything else in UFL."""
        return self.reuse_if_untouched(
            o,
            reference_value=reference_value,
            reference_grad=reference_grad,
            restricted=restricted,
        )

    @property
    def sub_graphs(self)-> list[ufl.core.expr.Expr]:
        return self._sub_graphs
