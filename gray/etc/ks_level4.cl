real gray_suba = a_spin*a_spin;
real gray_subb = z*z;
real gray_subc = gray_suba*gray_subb;
real gray_subd = x*x;
real gray_sube = y*y;
real gray_subf = gray_subb + gray_subd + gray_sube;
real gray_subg = -gray_suba + gray_subf;
real gray_subh = GRAY_SQRT(gray_subc + K(0.25)*(gray_subg*gray_subg));
real gray_subi = K(0.5)*gray_suba;
real gray_subj = gray_subh - gray_subi;
real gray_subba = K(0.5)*gray_subb + K(0.5)*gray_subd + K(0.5)*gray_sube;
real gray_subbb = gray_subba + gray_subj;
real gray_subbc = gray_subc + gray_subbb*gray_subbb;
real gray_subbd = K(1.0)/gray_subbc;
real gray_subbe = K(2)*x;
real gray_subbf = K(1.0)*x;
real gray_subbg = K(1.0)/gray_subh;
real gray_subbh = gray_subbg*gray_subg;
real gray_subbi = gray_subbf*gray_subbh;
real gray_subbj = gray_subbd*(gray_subbe + gray_subbi);
real gray_subca = gray_subbj*z;
real gray_subcb = K(2.0)*x;
real gray_subcc = gray_subbh*x;
real gray_subcd = gray_subcb + gray_subcc;
real gray_subce = K(-1.0)*gray_suba + gray_subf + K(2)*gray_subh;
real gray_subcf = K(1.0)/(gray_subbc*gray_subbc);
real gray_subcg = gray_subce*gray_subcf;
real gray_subch = gray_subcd*gray_subcg;
real gray_subci = gray_subbb*z;
real gray_subcj = gray_subch*gray_subci;
real gray_subda = gray_subca - gray_subcj;
real gray_subdb = GRAY_SQRT(gray_subbb);
real gray_subdc = gray_subbj*gray_subdb;
real gray_subdd = gray_subbd*gray_subce;
real gray_subde = K(1.0)/gray_subdb;
real gray_subdf = K(0.5)*x;
real gray_subdg = K(0.25)*gray_subcc;
real gray_subdh = gray_subde*(gray_subdf + gray_subdg);
real gray_subdi = gray_subdd*gray_subdh;
real gray_subdj = GRAY_SQRT_CUBE(gray_subbb);
real gray_subea = gray_subcg*gray_subdj;
real gray_subeb = gray_subcd*gray_subea;
real gray_subec = -a_spin*x + gray_subdb*y;
real gray_subed = gray_subba + gray_subh + gray_subi;
real gray_subee = K(1.0)/gray_subed;
real gray_subef = gray_subdc*gray_subee;
real gray_subeg = -a_spin + gray_subdh*y;
real gray_subeh = gray_subdb*gray_subdd;
real gray_subei = gray_subee*gray_subeh;
real gray_subej = -gray_subbf - gray_subbh*gray_subdf;
real gray_subfa = K(1.0)/(gray_subed*gray_subed);
real gray_subfb = gray_subeh*gray_subfa;
real gray_subfc = gray_subej*gray_subfb;
real gray_subfd = gray_subec*gray_subee;
real gray_subfe = gray_subdi*gray_subfd - gray_subeb*gray_subfd + gray_subec*gray_subef + gray_subec*gray_subfc + gray_subeg*gray_subei;
real gray_subff = a_spin*y + gray_subdb*x;
real gray_subfg = gray_subee*gray_subff;
real gray_subfh = gray_subdb + gray_subdh*x;
real gray_subfi = gray_subdi*gray_subfg - gray_subeb*gray_subfg + gray_subef*gray_subff + gray_subei*gray_subfh + gray_subfc*gray_subff;
real gray_subfj = gray_subca*gray_subee;
real gray_subga = gray_subdd*z;
real gray_subgb = gray_subfa*gray_subga;
real gray_subgc = gray_subej*gray_subgb;
real gray_subgd = gray_subdd*gray_subee;
real gray_subge = gray_subgd*z;
real gray_subgf = -gray_subcj*gray_subfg + gray_subff*gray_subfj + gray_subff*gray_subgc + gray_subfh*gray_subge;
real gray_subgg = gray_subff*gray_subff;
real gray_subgh = gray_subdc*gray_subfa;
real gray_subgi = K(1.0)/(gray_subed*gray_subed*gray_subed);
real gray_subgj = gray_subeh*gray_subgi;
real gray_subha = gray_subgj*(-gray_subbi - gray_subcb);
real gray_subhb = gray_subfa*gray_subgg;
real gray_subhc = K(2)*gray_subdb;
real gray_subhd = gray_subfb*gray_subff;
real gray_subhe = gray_subec*gray_subff;
real gray_subhf = gray_subfa*gray_subhe;
real gray_subhg = gray_subec*gray_subfb;
real gray_subhh = gray_subdi*gray_subhf - gray_subeb*gray_subhf + gray_subeg*gray_subhd + gray_subfh*gray_subhg + gray_subgh*gray_subhe + gray_subha*gray_subhe;
real gray_subhi = -gray_subcj*gray_subfd + gray_subec*gray_subfj + gray_subec*gray_subgc + gray_subeg*gray_subge;
real gray_subhj = gray_subec*gray_subec;
real gray_subia = gray_subfa*gray_subhj;
real gray_subib = K(2)*a_spin;
real gray_subic = K(2)*y;
real gray_subid = gray_subb*gray_subde;
real gray_subie = gray_subb*gray_subdd*K(1.0)/gray_subdj;
real gray_subif = gray_subb*gray_subdb;
real gray_subig = K(1.0)*y;
real gray_subih = gray_subbh*gray_subig;
real gray_subii = gray_subbd*(gray_subic + gray_subih);
real gray_subij = gray_subii*z;
real gray_subja = K(2.0)*y;
real gray_subjb = gray_subbh*y;
real gray_subjc = gray_subja + gray_subjb;
real gray_subjd = gray_subcg*gray_subjc;
real gray_subje = gray_subci*gray_subjd;
real gray_subjf = gray_subij - gray_subje;
real gray_subjg = gray_subdb*gray_subii;
real gray_subjh = K(0.5)*y;
real gray_subji = K(0.25)*gray_subjb;
real gray_subjj = gray_subde*(gray_subjh + gray_subji);
real gray_subbaa = gray_subdd*gray_subjj;
real gray_subbab = gray_subea*gray_subjc;
real gray_subbac = gray_subee*gray_subjg;
real gray_subbad = a_spin + gray_subjj*x;
real gray_subbae = -gray_subbh*gray_subjh - gray_subig;
real gray_subbaf = gray_subbaa*gray_subfg - gray_subbab*gray_subfg + gray_subbac*gray_subff + gray_subbad*gray_subei + gray_subbae*gray_subhd;
real gray_subbag = gray_subdb + gray_subjj*y;
real gray_subbah = gray_subbaa*gray_subfd - gray_subbab*gray_subfd + gray_subbac*gray_subec + gray_subbae*gray_subhg + gray_subbag*gray_subei;
real gray_subbai = gray_subee*gray_subij;
real gray_subbaj = gray_subbae*gray_subgb;
real gray_subbba = gray_subbad*gray_subge + gray_subbai*gray_subff + gray_subbaj*gray_subff - gray_subfg*gray_subje;
real gray_subbbb = gray_subfa*gray_subjg;
real gray_subbbc = gray_subgj*(-gray_subih - gray_subja);
real gray_subbbd = gray_subbaa*gray_subhf - gray_subbab*gray_subhf + gray_subbad*gray_subhg + gray_subbag*gray_subhd + gray_subbbb*gray_subhe + gray_subbbc*gray_subhe;
real gray_subbbe = gray_subbag*gray_subge + gray_subbai*gray_subec + gray_subbaj*gray_subec - gray_subfd*gray_subje;
real gray_subbbf = K(2)*z;
real gray_subbbg = K(0.5)*z;
real gray_subbbh = gray_subbg*(gray_suba*z + gray_subbbg*gray_subg);
real gray_subbbi = K(2)*gray_subbbh;
real gray_subbbj = gray_subbd*(gray_subbbf + gray_subbbi);
real gray_subbca = gray_subbbj*z;
real gray_subbcb = K(2.0)*z;
real gray_subbcc = gray_subcg*(-gray_suba*gray_subbbf - gray_subbb*(gray_subbbi + gray_subbcb));
real gray_subbcd = gray_subbcc*z;
real gray_subbce = gray_subbca + gray_subbcd + gray_subdd;
real gray_subbcf = gray_subbbj*gray_subdb;
real gray_subbcg = (K(1.0)/K(2.0))*gray_subbbh;
real gray_subbch = gray_subdd*(gray_subbbg + gray_subbcg);
real gray_subbci = gray_subbch*gray_subde;
real gray_subbcj = gray_subbcc*gray_subdb;
real gray_subbda = gray_subbch*gray_subee;
real gray_subbdb = gray_subbcf*gray_subee;
real gray_subbdc = -gray_subbbh - K(1.0)*z;
real gray_subbdd = gray_subbci*gray_subfg + gray_subbcj*gray_subfg + gray_subbda*x + gray_subbdb*gray_subff + gray_subbdc*gray_subhd;
real gray_subbde = gray_subbci*gray_subfd + gray_subbcj*gray_subfd + gray_subbda*y + gray_subbdb*gray_subec + gray_subbdc*gray_subhg;
real gray_subbdf = gray_subff*gray_subgd;
real gray_subbdg = gray_subbca*gray_subee;
real gray_subbdh = gray_subbci*gray_subee*z;
real gray_subbdi = gray_subbdc*gray_subgb;
real gray_subbdj = gray_subbcd*gray_subfg + gray_subbdf + gray_subbdg*gray_subff + gray_subbdh*x + gray_subbdi*gray_subff;
real gray_subbea = gray_subbch*gray_subfa;
real gray_subbeb = gray_subbea*gray_subff;
real gray_subbec = gray_subbcf*gray_subfa;
real gray_subbed = gray_subgj*(-gray_subbbi - gray_subbcb);
real gray_subbee = gray_subbea*gray_subec;
real gray_subbef = gray_subbci*gray_subhf + gray_subbcj*gray_subhf + gray_subbeb*y + gray_subbec*gray_subhe + gray_subbed*gray_subhe + gray_subbee*x;
real gray_subbeg = gray_subec*gray_subgd;
real gray_subbeh = gray_subbcd*gray_subfd + gray_subbdg*gray_subec + gray_subbdh*y + gray_subbdi*gray_subec + gray_subbeg;
real gray_subbei = K(1.0)/(gray_subed*gray_subed*gray_subed*gray_subed);
real gray_subbej = (K(1.0)/K(2.0))*gray_subb + (K(1.0)/K(2.0))*gray_subd + (K(1.0)/K(2.0))*gray_sube + gray_subj;
real gray_subbfa = K(8)*(gray_subbej*gray_subbej*gray_subbej);
real gray_subbfb = K(1.0)/(gray_subbc*gray_subbc*gray_subbc);
real gray_subbfc = K(4)*(gray_subbej*gray_subbej);
real gray_subbfd = gray_subdd*gray_subid + K(1);
real gray_subbfe = gray_subfb*gray_subgg + K(1);
real gray_subbff = gray_subfb*gray_subhj + K(1);
real gray_subbfg = gray_subbfd*gray_subbff;
real gray_subbfh = gray_subbfe*gray_subbfg;
real gray_subbfi = gray_subeh - K(1);
real gray_subbfj = K(1.0)/(K(-48)*gray_subb*gray_subbb*gray_subbei*gray_subgg*gray_subhj*K(1.0)/(gray_subbc*gray_subbc*gray_subbc*gray_subbc)*gray_subbej*gray_subbej*gray_subbej*gray_subbej + K(2)*gray_subb*gray_subbei*gray_subbfa*gray_subbfb*gray_subbfi*gray_subdb*gray_subgg*gray_subhj + K(2)*gray_subb*gray_subbfa*gray_subbfb*gray_subbfe*gray_subdb*gray_subfa*gray_subhj + K(2)*gray_subb*gray_subbfa*gray_subbfb*gray_subbff*gray_subdb*gray_subfa*gray_subgg - gray_subb*gray_subbfc*gray_subbfe*gray_subbff*gray_subcf - gray_subb*gray_subbfc*gray_subbfe*gray_subbfi*gray_subcf*gray_subfa*gray_subhj - gray_subb*gray_subbfc*gray_subbff*gray_subbfi*gray_subcf*gray_subfa*gray_subgg - gray_subbb*gray_subbei*gray_subbfc*gray_subbfd*gray_subbfi*gray_subcf*gray_subgg*gray_subhj - gray_subbb*gray_subbfc*gray_subbfd*gray_subbfe*gray_subcf*gray_subfa*gray_subhj - gray_subbb*gray_subbfc*gray_subbfd*gray_subbff*gray_subcf*gray_subfa*gray_subgg + K(2)*gray_subbei*gray_subbfa*gray_subbfb*gray_subbfd*gray_subdj*gray_subgg*gray_subhj + gray_subbfh*gray_subbfi);
real gray_subbga = gray_subbfj*(-gray_subb*gray_subbfa*gray_subbfb*gray_subdb*gray_subff*gray_subgi*gray_subhj + gray_subb*gray_subbfc*gray_subbff*gray_subcf*gray_subee*gray_subff + gray_subbb*gray_subbfc*gray_subbfd*gray_subcf*gray_subff*gray_subgi*gray_subhj - gray_subbfg*gray_subei*gray_subff);
real gray_subbgb = gray_subbfj*(-gray_subb*gray_subbfa*gray_subbfb*gray_subdb*gray_subec*gray_subgg*gray_subgi + gray_subb*gray_subbfc*gray_subbfe*gray_subcf*gray_subec*gray_subee + gray_subbb*gray_subbfc*gray_subbfd*gray_subcf*gray_subec*gray_subgg*gray_subgi - gray_subbfd*gray_subbfe*gray_subec*gray_subei);
real gray_subbgc = gray_subbfj*(-gray_subbb*gray_subbei*gray_subbfa*gray_subbfb*gray_subgg*gray_subhj*z + gray_subbfc*gray_subbfe*gray_subcf*gray_subdb*gray_subfa*gray_subhj*z + gray_subbfc*gray_subbff*gray_subcf*gray_subdb*gray_subfa*gray_subgg*z - gray_subbfe*gray_subbff*gray_subga);
real gray_subbgd = gray_subbfj*(-gray_subb*gray_subbfa*gray_subbfb*gray_subdb*gray_subec*gray_subfa*gray_subff + gray_subb*gray_subbfc*gray_subbfi*gray_subcf*gray_subec*gray_subfa*gray_subff + gray_subbb*gray_subbfc*gray_subbfd*gray_subcf*gray_subec*gray_subfa*gray_subff - gray_subbfd*gray_subbfi*gray_subec*gray_subhd);
real gray_subbge = gray_subbfj*(-gray_subbb*gray_subbfa*gray_subbfb*gray_subff*gray_subgi*gray_subhj*z - gray_subbdf*gray_subbff*gray_subbfi*z + gray_subbfc*gray_subbff*gray_subcf*gray_subdb*gray_subee*gray_subff*z + gray_subbfc*gray_subbfi*gray_subcf*gray_subdb*gray_subff*gray_subgi*gray_subhj*z);
real gray_subbgf = gray_subbfe*gray_subbfi;
real gray_subbgg = gray_subbfj*(-gray_subbb*gray_subbfa*gray_subbfb*gray_subec*gray_subgg*gray_subgi*z - gray_subbeg*gray_subbgf*z + gray_subbfc*gray_subbfe*gray_subcf*gray_subdb*gray_subec*gray_subee*z + gray_subbfc*gray_subbfi*gray_subcf*gray_subdb*gray_subec*gray_subgg*gray_subgi*z);
real16 Upsilon = {K(0),
K(0),
K(0),
K(0),
-gray_subda*uUPz - gray_subfe*uUPy - gray_subfi*uUPx - uUPt*(gray_subdc + gray_subdi - gray_subeb),
-gray_subfi*uUPt - gray_subgf*uUPz - gray_subhh*uUPy - uUPx*(gray_subdi*gray_subhb - gray_subeb*gray_subhb + gray_subgg*gray_subgh + gray_subgg*gray_subha + gray_subhd*(gray_subbe*gray_subdh + gray_subhc)),
-gray_subfe*uUPt - gray_subhh*uUPx - gray_subhi*uUPz - uUPy*(gray_subdi*gray_subia - gray_subeb*gray_subia + gray_subgh*gray_subhj + gray_subha*gray_subhj + gray_subhg*(gray_subdh*gray_subic - gray_subib)),
-gray_subda*uUPt - gray_subgf*uUPx - gray_subhi*uUPy - uUPz*(gray_subbj*gray_subid - gray_subch*gray_subif + gray_subie*(-gray_subdf - gray_subdg)),
-gray_subbaf*uUPx - gray_subbah*uUPy - gray_subjf*uUPz - uUPt*(gray_subbaa - gray_subbab + gray_subjg),
-gray_subbaf*uUPt - gray_subbba*uUPz - gray_subbbd*uUPy - uUPx*(gray_subbaa*gray_subhb - gray_subbab*gray_subhb + gray_subbbb*gray_subgg + gray_subbbc*gray_subgg + gray_subhd*(gray_subbe*gray_subjj + gray_subib)),
-gray_subbah*uUPt - gray_subbbd*uUPx - gray_subbbe*uUPz - uUPy*(gray_subbaa*gray_subia - gray_subbab*gray_subia + gray_subbbb*gray_subhj + gray_subbbc*gray_subhj + gray_subhg*(gray_subhc + gray_subic*gray_subjj)),
-gray_subbba*uUPx - gray_subbbe*uUPy - gray_subjf*uUPt - uUPz*(gray_subid*gray_subii + gray_subie*(-gray_subjh - gray_subji) - gray_subif*gray_subjd),
-gray_subbce*uUPz - gray_subbdd*uUPx - gray_subbde*uUPy - uUPt*(gray_subbcf + gray_subbci + gray_subbcj),
-gray_subbdd*uUPt - gray_subbdj*uUPz - gray_subbef*uUPy - uUPx*(gray_subbci*gray_subhb + gray_subbcj*gray_subhb + gray_subbe*gray_subbeb + gray_subbec*gray_subgg + gray_subbed*gray_subgg),
-gray_subbde*uUPt - gray_subbef*uUPx - gray_subbeh*uUPz - uUPy*(gray_subbci*gray_subia + gray_subbcj*gray_subia + gray_subbec*gray_subhj + gray_subbed*gray_subhj + gray_subbee*gray_subic),
-gray_subbce*uUPt - gray_subbdj*uUPx - gray_subbeh*uUPy - uUPz*(gray_subbbf*gray_subdd*gray_subde + gray_subbbj*gray_subid + gray_subbcc*gray_subid + gray_subie*(-gray_subbbg - gray_subbcg))};
real16 gUP = {gray_subbfj*(K(2)*gray_subb*gray_subbei*gray_subbfa*gray_subbfb*gray_subdb*gray_subgg*gray_subhj - gray_subb*gray_subbfc*gray_subbfe*gray_subcf*gray_subfa*gray_subhj - gray_subb*gray_subbfc*gray_subbff*gray_subcf*gray_subfa*gray_subgg - gray_subbb*gray_subbei*gray_subbfc*gray_subbfd*gray_subcf*gray_subgg*gray_subhj + gray_subbfh),
gray_subbga,
gray_subbgb,
gray_subbgc,
gray_subbga,
gray_subbfj*(K(2)*gray_subb*gray_subbfa*gray_subbfb*gray_subdb*gray_subfa*gray_subhj - gray_subb*gray_subbfc*gray_subbff*gray_subcf - gray_subb*gray_subbfc*gray_subbfi*gray_subcf*gray_subfa*gray_subhj - gray_subbb*gray_subbfc*gray_subbfd*gray_subcf*gray_subfa*gray_subhj + gray_subbfg*gray_subbfi),
gray_subbgd,
gray_subbge,
gray_subbgb,
gray_subbgd,
gray_subbfj*(K(2)*gray_subb*gray_subbfa*gray_subbfb*gray_subdb*gray_subfa*gray_subgg - gray_subb*gray_subbfc*gray_subbfe*gray_subcf - gray_subb*gray_subbfc*gray_subbfi*gray_subcf*gray_subfa*gray_subgg - gray_subbb*gray_subbfc*gray_subbfd*gray_subcf*gray_subfa*gray_subgg + gray_subbfd*gray_subbgf),
gray_subbgg,
gray_subbgc,
gray_subbge,
gray_subbgg,
gray_subbfj*(-gray_subbb*gray_subbei*gray_subbfc*gray_subbfi*gray_subcf*gray_subgg*gray_subhj - gray_subbb*gray_subbfc*gray_subbfe*gray_subcf*gray_subfa*gray_subhj - gray_subbb*gray_subbfc*gray_subbff*gray_subcf*gray_subfa*gray_subgg + K(2)*gray_subbei*gray_subbfa*gray_subbfb*gray_subdj*gray_subgg*gray_subhj + gray_subbff*gray_subbgf)};
real16 UpsilonT = Upsilon.s048c159d26ae37bf;
real16 Xi = UpsilonT - Upsilon/2;
real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));
