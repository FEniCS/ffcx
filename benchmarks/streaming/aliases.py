"""Experimental propagation of scalar aliases into contraction grouping."""

# Generated C fragments retain their original source layout.

import ffcx.codegeneration.lnodes as L


def install():
    """Install the experimental generator monkeypatch."""
    from ffcx.codegeneration.integral_generator import IntegralGenerator

    original_partition = IntegralGenerator.generate_partition
    original_get = IntegralGenerator.get_var

    def get(self, *args, **kwargs):
        """Apply one step of the experimental transformation."""
        result = original_get(self, *args, **kwargs)
        aliases = getattr(self, "_aliases", {})
        while isinstance(result, L.Symbol) and result.name in aliases:
            result = aliases[result.name]
        return result

    def partition(self, *args, **kwargs):
        """Apply one step of the experimental transformation."""
        definitions, intermediates = original_partition(self, *args, **kwargs)
        if not hasattr(self, "_aliases"):
            self._aliases = {}
        for s in intermediates:
            if (
                isinstance(s, L.VariableDecl)
                and isinstance(s.value, L.Symbol)
                and s.symbol.dtype == s.value.dtype
            ):
                value = s.value
                while value.name in self._aliases:
                    value = self._aliases[value.name]
                self._aliases[s.symbol.name] = value
        return definitions, intermediates

    IntegralGenerator.get_var = get
    IntegralGenerator.generate_partition = partition
