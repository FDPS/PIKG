def add_build_in_functions(main_code):
    code = ""
    code += rsqrt_fnc
    code += inv_fnc
    code += main_code # main
    return code

#1次精度を2iter.
rsqrt_fnc = """
function rsqrt(x)
    y0 = rsqrt_approx(x)
    y1 = y0 + y0*(0.5-0.5*(y0*y0*x))
    y2 = y1 + y1*(0.5-0.5*(y1*y1*x))
    return y2
end
"""

#1次精度を2iter.
inv_fnc = """
function inv(x)
    y0r = rsqrt_approx(x)
    y0 = y0r * y0r
    y1 = y0*(2.0-(y0*x))
    y2 = y1*(2.0-(y1*x))
    return y2
end
"""
    
    
