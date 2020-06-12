real f,  dx_f,  dy_f,  dz_f;
real lx, dx_lx, dy_lx, dz_lx;
real ly, dx_ly, dy_ly, dz_ly;
real lz, dx_lz, dy_lz, dz_lz;

real  hDxu, hDyu, hDzu;
real4 uD;
real  tmp;

{
    real dx_r, dy_r, dz_r;
    real r, ir, iss;
    {
        real aa = a_spin * a_spin;
        real rr;
        {
            real zz = q.s3 * q.s3;
            real dd;
            {
                real kk = K(0.5) * (q.s1 * q.s1 + q.s2 * q.s2 + zz - aa);
                dd = sqrt(kk * kk + aa * zz);
                rr = dd + kk;
            }
            r  = sqrt(rr);
            printf("\n r = %.16g ", r);
            ir = K(1.0) / r;
            {
                real ss = rr + aa;
                iss  = K(1.0) / ss;
                tmp  = K(0.5) / (r * dd);
                dz_r = tmp * ss * q.s3;
                tmp *= rr;
            }
            dy_r = tmp * q.s2;
            dx_r = tmp * q.s1;
            tmp  = K(2.0) / (rr + aa * zz / rr);
        }
        f    = tmp *  r;
        dx_f = tmp *  dx_r * (K(3.0) - K(2.0) * rr * tmp);
        dy_f = tmp *  dy_r * (K(3.0) - K(2.0) * rr * tmp);
        dz_f = tmp * (dz_r * (K(3.0) - K(2.0) * rr * tmp) - tmp * aa * q.s3 * ir);
    } /* 48 (-8) FLOPs; estimated FLoating-point OPerations, the number
         in the parentheses is (the negative of) the number of FMA */
    {
        real m2r  = K(-2.0) * r;
        real issr =     iss * r;
        real issa =     iss * a_spin;

        lx    = iss * (q.s1 * r + q.s2 * a_spin);
        tmp   = iss * (q.s1 + m2r * lx);
        dx_lx = tmp * dx_r + issr;
        dy_lx = tmp * dy_r + issa;
        dz_lx = tmp * dz_r;

        ly    = iss * (q.s2 * r - q.s1 * a_spin);
        tmp   = iss * (q.s2 + m2r * ly);
        dx_ly = tmp * dx_r - issa;
        dy_ly = tmp * dy_r + issr;
        dz_ly = tmp * dz_r;

        lz    = q.s3 * ir;
        tmp   = -lz * ir;
        dx_lz = tmp * dx_r;
        dy_lz = tmp * dy_r;
        dz_lz = tmp * dz_r + ir;
    } /* 35 (-9) FLOPs */
}

{
    real  flu;
    real4 Dx, Dy, Dz;
    {
        real lu = u.s0 + lx * u.s1 + ly * u.s2 + lz * u.s3;
        flu   = f * lu;
        Dx.s0 = dx_f * lu + f * (dx_lx * u.s1 + dx_ly * u.s2 + dx_lz * u.s3);
        Dy.s0 = dy_f * lu + f * (dy_lx * u.s1 + dy_ly * u.s2 + dy_lz * u.s3);
        Dz.s0 = dz_f * lu + f * (dz_lx * u.s1 + dz_ly * u.s2 + dz_lz * u.s3);
    }
    Dx.s1 = Dx.s0 * lx + flu * dx_lx;
    Dx.s2 = Dx.s0 * ly + flu * dx_ly;
    Dx.s3 = Dx.s0 * lz + flu * dx_lz; /* 9 (-3) FLOPs */

    Dy.s1 = Dy.s0 * lx + flu * dy_lx;
    Dy.s2 = Dy.s0 * ly + flu * dy_ly;
    Dy.s3 = Dy.s0 * lz + flu * dy_lz; /* 9 (-3) FLOPs */

    Dz.s1 = Dz.s0 * lx + flu * dz_lx;
    Dz.s2 = Dz.s0 * ly + flu * dz_ly;
    Dz.s3 = Dz.s0 * lz + flu * dz_lz; /* 9 (-3) FLOPs */

    hDxu = K(0.5) * dot(Dx, u);
    hDyu = K(0.5) * dot(Dy, u);
    hDzu = K(0.5) * dot(Dz, u); /* 24 (-9) FLOPs */

    uD  = u.s1 * Dx + u.s2 * Dy + u.s3 * Dz; /* 20 (-8) FLOPs */

    tmp = f * (-uD.s0 + lx * (uD.s1 - hDxu) + ly * (uD.s2 - hDyu) + lz * (uD.s3 - hDzu));
}
printf("\n DIFF 0 %.16g %.16g %.16g ",
        fabs((rhs.s0 - (uD.s0 - tmp))/(uD.s0 - tmp)),
        rhs.s0, (uD.s0 - tmp));
printf("\n DIFF 1 %.16g %.16g %.16g ",
        fabs((rhs.s1 - (hDxu - uD.s1 + lx * tmp))/(hDxu - uD.s1 + lx * tmp)),
        rhs.s1, (hDxu - uD.s1 + lx * tmp));
printf("\n DIFF 2 %.16g %.16g %.16g ",
        fabs((rhs.s2 - (hDyu - uD.s2 + ly * tmp))/(hDyu - uD.s2 + ly * tmp)),
        rhs.s2, (hDyu - uD.s2 + ly * tmp));
printf("\n DIFF 3 %.16g %.16g %.16g ",
        fabs((rhs.s3 - (hDzu - uD.s3 + lz * tmp))/(hDzu - uD.s3 + lz * tmp)),
        rhs.s3, (hDzu - uD.s3 + lz * tmp));
