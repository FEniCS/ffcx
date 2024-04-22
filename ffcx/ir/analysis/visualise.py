# Copyright (C) 2018 Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utility to draw graphs."""

from ufl.classes import (
    Argument,
    Division,
    FloatValue,
    Indexed,
    IntValue,
    Product,
    ReferenceValue,
    Sum,
)

from ffcx.ir.analysis.modified_terminals import strip_modified_terminal


def visualise_graph(Gx, filename):
    """Visualise a graph."""
    try:
        import pygraphviz as pgv
    except ImportError:
        raise RuntimeError("Install pygraphviz")

    if Gx.number_of_nodes() > 400:
        print("Skipping visualisation")
        return

    G = pgv.AGraph(strict=False, directed=True)
    for nd, v in Gx.nodes.items():
        ex = v["expression"]
        label = ex.__class__.__name__
        if isinstance(ex, Sum):
            label = "+"
        elif isinstance(ex, Product):
            label = "*"
        elif isinstance(ex, Division):
            label = "/"
        elif isinstance(ex, (IntValue, FloatValue)):
            label = ex.value()
        elif isinstance(ex, (Indexed, ReferenceValue)):
            label = str(ex)
        G.add_node(nd, label="[%d] %s" % (nd, label))

        arg = strip_modified_terminal(ex)
        if isinstance(arg, Argument):
            G.get_node(nd).attr["shape"] = "box"

        stat = v.get("status")
        if stat == "piecewise":
            G.get_node(nd).attr["color"] = "blue"
            G.get_node(nd).attr["penwidth"] = 5
        elif stat == "varying":
            G.get_node(nd).attr["color"] = "red"
            G.get_node(nd).attr["penwidth"] = 5
        elif stat == "inactive":
            G.get_node(nd).attr["color"] = "dimgray"
            G.get_node(nd).attr["penwidth"] = 5

        t = v.get("target")
        if t:
            G.get_node(nd).attr["label"] += ":" + str(t)
            G.get_node(nd).attr["shape"] = "hexagon"

        c = v.get("component")
        if c:
            G.get_node(nd).attr["label"] += f", comp={c}"

    for nd, eds in Gx.out_edges.items():
        for ed in eds:
            G.add_edge(nd, ed)

    G.layout(prog="dot")
    G.draw(filename)
