import sympy as sp

def verify_C_implementation():
    print("=" * 65)
    print(" SYMBOLIC VERIFICATION OF ZONAL HARMONICS (J6, J8, J10)")
    print("=" * 65)

    x, y, z = sp.symbols('x y z', real=True)
    G, m, R_eq = sp.symbols('G m R_eq', real=True, positive=True)
    J6, J8, J10 = sp.symbols('J6 J8 J10', real=True)

    r = sp.sqrt(x**2 + y**2 + z**2)
    u = z / r
    costheta2 = u**2

    armonicos = [
        (6, J6, 
         sp.Rational(7, 16) * G * m * J6 * (R_eq**6) / (r**9),
         429 * (costheta2**3) - 495 * (costheta2**2) + 135 * costheta2 - 5,
         429 * (costheta2**3) - 693 * (costheta2**2) + 315 * costheta2 - 35),
        
        (8, J8, 
         sp.Rational(9, 128) * G * m * J8 * (R_eq**8) / (r**11),
         12155 * (costheta2**4) - 20020 * (costheta2**3) + 10010 * (costheta2**2) - 1540 * costheta2 + 35,
         12155 * (costheta2**4) - 25740 * (costheta2**3) + 18018 * (costheta2**2) - 4620 * costheta2 + 315),
        
        (10, J10, 
         sp.Rational(11, 256) * G * m * J10 * (R_eq**10) / (r**13),
         88179 * (costheta2**5) - 188955 * (costheta2**4) + 139230 * (costheta2**3) - 40950 * (costheta2**2) + 4095 * costheta2 - 63,
         88179 * (costheta2**5) - 230945 * (costheta2**4) + 218790 * (costheta2**3) - 90090 * (costheta2**2) + 15015 * costheta2 - 693)
    ]

    for n, J_sym, f1_code, f2_code, f3_code in armonicos:
        print(f"\n[+] Harmonic evaluation J{n} (calculating the gradient 3D...).")
        
        Pn = sp.legendre(n, u)
        V_pert = (G * m * (R_eq**n) * J_sym / (r**(n+1))) * Pn
        
        ax_exact = -sp.diff(V_pert, x)
        ay_exact = -sp.diff(V_pert, y)
        az_exact = -sp.diff(V_pert, z)
        
        ax_code = f1_code * f2_code * x
        ay_code = f1_code * f2_code * y
        az_code = f1_code * f3_code * z
        
        diff_x = sp.simplify(ax_exact - ax_code)
        diff_y = sp.simplify(ay_exact - ay_code)
        diff_z = sp.simplify(az_exact - az_code)

        assert diff_x == 0
        assert diff_y == 0
        assert diff_z == 0

verify_C_implementation()
