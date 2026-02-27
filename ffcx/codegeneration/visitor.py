# Copyright (C) 2026
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Visitor and Transformer base classes for LNode AST traversal."""

from ffcx.codegeneration.lnodes import LNode


class LNodeVisitor:
    """Read-only visitor for LNode trees.

    Dispatches to ``visit_ClassName`` methods based on the node type.
    Override specific ``visit_*`` methods to handle particular node types.
    The default ``generic_visit`` recurses into all children.
    """

    def visit(self, node: LNode):
        """Dispatch to the appropriate visit method for *node*."""
        method_name = f"visit_{type(node).__name__}"
        visitor = getattr(self, method_name, self.generic_visit)
        return visitor(node)

    def generic_visit(self, node: LNode):
        """Visit all children of *node*."""
        for child in node.children():
            self.visit(child)


class LNodeTransformer:
    """Functional transformer for LNode trees.

    Like ``LNodeVisitor``, but each ``visit_*`` method should return an LNode.
    The default ``generic_visit`` rebuilds the node with transformed children.
    If no children change, the original node is returned (no unnecessary copies).
    """

    def visit(self, node: LNode) -> LNode:
        """Dispatch to the appropriate visit method for *node*."""
        method_name = f"visit_{type(node).__name__}"
        transformer = getattr(self, method_name, self.generic_visit)
        return transformer(node)

    def generic_visit(self, node: LNode) -> LNode:
        """Rebuild *node* with transformed children.

        Returns the original node if no children changed.
        """
        named = node.named_children()
        if not named:
            return node

        changes: dict[str, object] = {}
        for field_name, child in named:
            if isinstance(child, list):
                new_list = [self.visit(c) if isinstance(c, LNode) else c for c in child]
                if any(a is not b for a, b in zip(child, new_list)):
                    changes[field_name] = new_list
            elif isinstance(child, tuple):
                new_tuple = tuple(self.visit(c) if isinstance(c, LNode) else c for c in child)
                if any(a is not b for a, b in zip(child, new_tuple)):
                    changes[field_name] = new_tuple
            elif isinstance(child, LNode):
                new_child = self.visit(child)
                if new_child is not child:
                    changes[field_name] = new_child

        if not changes:
            return node
        return node.replace(**changes)
