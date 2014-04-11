from uflacs.codeutils.format_code import (format_code,
                                                    strip_trailing_whitespace,
                                                    WithKeywords, Class, Indented)

eval_template = """
virtual void eval(Array<double>& values,
                  const Array<double>& x) const
{
%(eval_body)s
}"""

eval_cell_template = """
virtual void eval(Array<double>& values,
                  const Array<double>& x,
                  const ufc::cell& cell) const
{
%(eval_body)s
}"""

c_template = "double %s;"
mf_template = "boost::shared_ptr<MeshFunction<double> > %s;"
gf_template = "boost::shared_ptr<GenericFunction> %s;"
f_template = "boost::shared_ptr<Function> %s;"

def format_dolfin_expression(classname="MyExpression",
                             shape=(),
                             eval_body="",
                             constants=(),
                             mesh_functions=(),
                             functions=(),
                             generic_functions=()):
    assert all(isinstance(x,str) for x in constants)
    assert all(isinstance(x,str) for x in mesh_functions)
    assert all(isinstance(x,str) for x in generic_functions)
    assert all(isinstance(x,str) for x in functions)
    # ...
    r = len(shape)
    if r == 0:
        constructors = []
    elif r == 1:
        constructors = ["%s(): Expression(%d) {}" % (classname, shape[0])]
    elif r == 2:
        constructors = ["%s(): Expression(%d, %d) {}" % (classname, shape[0], shape[1])]
    else:
        error
    # Make class body
    classbody = []
    classbody += constructors
    classbody += [WithKeywords(eval_template, { 'eval_body': Indented(eval_body) })]
    members = []
    members += [c_template  % f for f in constants]
    members += [mf_template % f for f in mesh_functions]
    members += [gf_template % f for f in generic_functions]
    members += [f_template  % f for f in functions]
    if members:
        classbody += ['']
        classbody += members
    code = Class(name=classname,
                 superclass="Expression",
                 public_body=classbody)
    return strip_trailing_whitespace(format_code(code))
