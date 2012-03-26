
def tabulate_tensor_prototype():

    def build_body01():
        body = []
        body += ["decl"]
        body += ["for ()"]
        decl = ["decl01;"]
        body += [Block(decl)]
        return body

    def build_loop():
        preloop = ["preloop"]
        postloop = ["postloop"]
        loopheader = "for ()"
        body = ["body;"]
        return [preloop, loopheader, Block(body), postloop]

    def build_body1():
        body = []
        body += ["decl"]
        body += ["for (i1)"]
        decl = ["decl1;"]
        body += [Block(decl)]
        return body

    def build_body0():
        body = []
        body += ["decl"]
        body += ["for (i1)"]
        decl = [build_body1()]
        body += [Block(decl)]
        return body

    def build_bodyx():
        body = []
        body += ["compute x from iq;"]
        body += ["compute quantities depending on x"]
        body += build_i0loop()
        body += build_i1loop()
        return body

    def build_i0loop():
        header = "for (int i0 = 0; i0 < n0; ++i0)"
        decl = build_body0()
        return [header, Block(decl)]

    def build_i1loop():
        header = "for (int i1 = 0; i1 < n1; ++i1)"
        decl = build_body1()
        return [header, Block(decl)]

    def build_body():
        body = []
        body += ["compute quantities depending on cell and coeffs but not x;"]
        body += ["for (int iq = 0; iq < nq; ++iq)"]
        decl = build_bodyx()
        body += [Block(decl)]
        return body

    print format_code(build_body())
