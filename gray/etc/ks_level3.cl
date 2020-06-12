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
real gray_subdb = gray_subda*uUPz;
real gray_subdc = GRAY_SQRT(gray_subbb);
real gray_subdd = gray_subbj*gray_subdc;
real gray_subde = gray_subbd*gray_subce;
real gray_subdf = K(1.0)/gray_subdc;
real gray_subdg = K(0.5)*x;
real gray_subdh = K(0.25)*gray_subcc;
real gray_subdi = gray_subdf*(gray_subdg + gray_subdh);
real gray_subdj = gray_subde*gray_subdi;
real gray_subea = GRAY_SQRT_CUBE(gray_subbb);
real gray_subeb = gray_subcg*gray_subea;
real gray_subec = gray_subcd*gray_subeb;
real gray_subed = uUPt*(gray_subdd + gray_subdj - gray_subec);
real gray_subee = -a_spin*x + gray_subdc*y;
real gray_subef = gray_subba + gray_subh + gray_subi;
real gray_subeg = K(1.0)/gray_subef;
real gray_subeh = gray_subdd*gray_subeg;
real gray_subei = -a_spin + gray_subdi*y;
real gray_subej = gray_subdc*gray_subde;
real gray_subfa = gray_subeg*gray_subej;
real gray_subfb = -gray_subbf - gray_subbh*gray_subdg;
real gray_subfc = K(1.0)/(gray_subef*gray_subef);
real gray_subfd = gray_subej*gray_subfc;
real gray_subfe = gray_subfb*gray_subfd;
real gray_subff = gray_subee*gray_subeg;
real gray_subfg = gray_subdj*gray_subff - gray_subec*gray_subff + gray_subee*gray_subeh + gray_subee*gray_subfe + gray_subei*gray_subfa;
real gray_subfh = gray_subfg*uUPy;
real gray_subfi = a_spin*y + gray_subdc*x;
real gray_subfj = gray_subeg*gray_subfi;
real gray_subga = gray_subdc + gray_subdi*x;
real gray_subgb = gray_subdj*gray_subfj - gray_subec*gray_subfj + gray_subeh*gray_subfi + gray_subfa*gray_subga + gray_subfe*gray_subfi;
real gray_subgc = gray_subgb*uUPx;
real gray_subgd = K(2)*y;
real gray_subge = K(1.0)*y;
real gray_subgf = gray_subbh*gray_subge;
real gray_subgg = gray_subbd*(gray_subgd + gray_subgf);
real gray_subgh = gray_subgg*z;
real gray_subgi = K(2.0)*y;
real gray_subgj = gray_subbh*y;
real gray_subha = gray_subgi + gray_subgj;
real gray_subhb = gray_subcg*gray_subha;
real gray_subhc = gray_subci*gray_subhb;
real gray_subhd = gray_subgh - gray_subhc;
real gray_subhe = gray_subhd*uUPz;
real gray_subhf = gray_subdc*gray_subgg;
real gray_subhg = K(0.5)*y;
real gray_subhh = K(0.25)*gray_subgj;
real gray_subhi = gray_subdf*(gray_subhg + gray_subhh);
real gray_subhj = gray_subde*gray_subhi;
real gray_subia = gray_subeb*gray_subha;
real gray_subib = uUPt*(gray_subhf + gray_subhj - gray_subia);
real gray_subic = gray_subeg*gray_subhf;
real gray_subid = a_spin + gray_subhi*x;
real gray_subie = -gray_subbh*gray_subhg - gray_subge;
real gray_subif = gray_subfd*gray_subie;
real gray_subig = gray_subfa*gray_subid + gray_subfi*gray_subic + gray_subfi*gray_subif + gray_subfj*gray_subhj - gray_subfj*gray_subia;
real gray_subih = gray_subig*uUPx;
real gray_subii = gray_subdc + gray_subhi*y;
real gray_subij = gray_subee*gray_subic + gray_subee*gray_subif + gray_subfa*gray_subii + gray_subff*gray_subhj - gray_subff*gray_subia;
real gray_subja = gray_subij*uUPy;
real gray_subjb = K(2)*z;
real gray_subjc = K(0.5)*z;
real gray_subjd = gray_subbg*(gray_suba*z + gray_subg*gray_subjc);
real gray_subje = K(2)*gray_subjd;
real gray_subjf = gray_subbd*(gray_subjb + gray_subje);
real gray_subjg = gray_subjf*z;
real gray_subjh = K(2.0)*z;
real gray_subji = gray_subcg*(-gray_suba*gray_subjb - gray_subbb*(gray_subje + gray_subjh));
real gray_subjj = gray_subji*z;
real gray_subbaa = gray_subde + gray_subjg + gray_subjj;
real gray_subbab = gray_subbaa*uUPz;
real gray_subbac = gray_subdc*gray_subjf;
real gray_subbad = (K(1.0)/K(2.0))*gray_subjd;
real gray_subbae = gray_subbad + gray_subjc;
real gray_subbaf = gray_subde*gray_subdf;
real gray_subbag = gray_subbae*gray_subbaf;
real gray_subbah = gray_subdc*gray_subji;
real gray_subbai = uUPt*(gray_subbac + gray_subbag + gray_subbah);
real gray_subbaj = gray_subde*gray_subeg;
real gray_subbba = gray_subbae*gray_subbaj;
real gray_subbbb = gray_subbac*gray_subeg;
real gray_subbbc = -gray_subjd - K(1.0)*z;
real gray_subbbd = gray_subbbc*gray_subfd;
real gray_subbbe = gray_subbag*gray_subfj + gray_subbah*gray_subfj + gray_subbba*x + gray_subbbb*gray_subfi + gray_subbbd*gray_subfi;
real gray_subbbf = gray_subbbe*uUPx;
real gray_subbbg = gray_subbag*gray_subff + gray_subbah*gray_subff + gray_subbba*y + gray_subbbb*gray_subee + gray_subbbd*gray_subee;
real gray_subbbh = gray_subbbg*uUPy;
real gray_subbbi = gray_subca*gray_subeg;
real gray_subbbj = gray_subfc*gray_subfi;
real gray_subbca = gray_subde*z;
real gray_subbcb = gray_subbca*gray_subfb;
real gray_subbcc = gray_subbaj*z;
real gray_subbcd = gray_subbbi*gray_subfi + gray_subbbj*gray_subbcb + gray_subbcc*gray_subga - gray_subcj*gray_subfj;
real gray_subbce = (K(1.0)/K(2.0))*uUPz;
real gray_subbcf = (K(1.0)/K(2.0))*uUPt;
real gray_subbcg = gray_subfi*gray_subfi;
real gray_subbch = gray_subdd*gray_subfc;
real gray_subbci = K(1.0)/(gray_subef*gray_subef*gray_subef);
real gray_subbcj = gray_subbci*gray_subej;
real gray_subbda = gray_subbcj*(-gray_subbi - gray_subcb);
real gray_subbdb = gray_subbcg*gray_subfc;
real gray_subbdc = K(2)*gray_subdc;
real gray_subbdd = gray_subfd*gray_subfi;
real gray_subbde = (K(1.0)/K(2.0))*uUPx;
real gray_subbdf = gray_subee*gray_subfi;
real gray_subbdg = gray_subbdf*gray_subfc;
real gray_subbdh = gray_subee*gray_subfd;
real gray_subbdi = gray_subbch*gray_subbdf + gray_subbda*gray_subbdf + gray_subbdd*gray_subei + gray_subbdg*gray_subdj - gray_subbdg*gray_subec + gray_subbdh*gray_subga;
real gray_subbdj = (K(1.0)/K(2.0))*uUPy;
real gray_subbea = gray_subbbj*gray_subbca;
real gray_subbeb = gray_subbcc*gray_subid + gray_subbea*gray_subie + gray_subfj*gray_subgh - gray_subfj*gray_subhc;
real gray_subbec = gray_subbeb*uUPz;
real gray_subbed = gray_subee*gray_subfc;
real gray_subbee = gray_subbbi*gray_subee + gray_subbcb*gray_subbed + gray_subbcc*gray_subei - gray_subcj*gray_subff;
real gray_subbef = gray_subbee*uUPz;
real gray_subbeg = gray_subig*uUPt;
real gray_subbeh = gray_subfg*uUPt;
real gray_subbei = -gray_subgf - gray_subgi;
real gray_subbej = gray_subbcg*gray_subbcj;
real gray_subbfa = K(2)*a_spin;
real gray_subbfb = uUPx*(gray_subbdb*gray_subhf + gray_subbdb*gray_subhj - gray_subbdb*gray_subia + gray_subbdd*(gray_subbe*gray_subhi + gray_subbfa) + gray_subbei*gray_subbej);
real gray_subbfc = gray_subee*gray_subee;
real gray_subbfd = gray_subbfc*gray_subfc;
real gray_subbfe = uUPy*(gray_subbch*gray_subbfc + gray_subbda*gray_subbfc + gray_subbdh*(-gray_subbfa + gray_subdi*gray_subgd) + gray_subbfd*gray_subdj - gray_subbfd*gray_subec);
real gray_subbff = gray_subbcj*gray_subbdf;
real gray_subbfg = gray_subbdd*gray_subii + gray_subbdg*gray_subhf + gray_subbdg*gray_subhj - gray_subbdg*gray_subia + gray_subbdh*gray_subid + gray_subbei*gray_subbff;
real gray_subbfh = gray_subbfg*uUPy;
real gray_subbfi = gray_subbdi*uUPx;
real gray_subbfj = gray_subda*uUPt;
real gray_subbga = gray_subb*gray_subdf;
real gray_subbgb = gray_subb*gray_subde*K(1.0)/gray_subea;
real gray_subbgc = gray_subb*gray_subdc;
real gray_subbgd = uUPz*(gray_subbga*gray_subbj + gray_subbgb*(-gray_subdg - gray_subdh) - gray_subbgc*gray_subch);
real gray_subbge = gray_subbee*uUPy;
real gray_subbgf = gray_subbcd*uUPx;
real gray_subbgg = gray_subbaj*gray_subfi;
real gray_subbgh = gray_subbag*gray_subeg*z;
real gray_subbgi = gray_subbbc*gray_subbea + gray_subbgg + gray_subbgh*x + gray_subfj*gray_subjg + gray_subfj*gray_subjj;
real gray_subbgj = gray_subbgi*uUPz;
real gray_subbha = gray_subbbe*uUPt;
real gray_subbhb = gray_subbae*gray_subde;
real gray_subbhc = gray_subbbj*gray_subbhb;
real gray_subbhd = -gray_subje - gray_subjh;
real gray_subbhe = uUPx*(gray_subbac*gray_subbdb + gray_subbag*gray_subbdb + gray_subbah*gray_subbdb + gray_subbe*gray_subbhc + gray_subbej*gray_subbhd);
real gray_subbhf = gray_subbed*gray_subbhb;
real gray_subbhg = gray_subbac*gray_subbdg + gray_subbag*gray_subbdg + gray_subbah*gray_subbdg + gray_subbff*gray_subbhd + gray_subbhc*y + gray_subbhf*x;
real gray_subbhh = gray_subbhg*uUPy;
real gray_subbhi = gray_subbca*gray_subbed;
real gray_subbhj = gray_subbcc*gray_subii + gray_subbhi*gray_subie + gray_subff*gray_subgh - gray_subff*gray_subhc;
real gray_subbia = gray_subbcj*gray_subbfc;
real gray_subbib = gray_subhd*uUPt;
real gray_subbic = uUPz*(gray_subbga*gray_subgg + gray_subbgb*(-gray_subhg - gray_subhh) - gray_subbgc*gray_subhb);
real gray_subbid = gray_subbeb*uUPx;
real gray_subbie = gray_subbhj*uUPy;
real gray_subbif = gray_subbaj*gray_subee;
real gray_subbig = gray_subbbc*gray_subbhi + gray_subbgh*y + gray_subbif + gray_subff*gray_subjg + gray_subff*gray_subjj;
real gray_subbih = gray_subbig*uUPz;
real gray_subbii = gray_subbbg*uUPt;
real gray_subbij = uUPy*(gray_subbac*gray_subbfd + gray_subbag*gray_subbfd + gray_subbah*gray_subbfd + gray_subbhd*gray_subbia + gray_subbhf*gray_subgd);
real gray_subbja = gray_subbhg*uUPx;
real gray_subbjb = K(1.0)/(gray_subef*gray_subef*gray_subef*gray_subef);
real gray_subbjc = (K(1.0)/K(2.0))*gray_subb + (K(1.0)/K(2.0))*gray_subd + (K(1.0)/K(2.0))*gray_sube + gray_subj;
real gray_subbjd = K(8)*(gray_subbjc*gray_subbjc*gray_subbjc);
real gray_subbje = K(1.0)/(gray_subbc*gray_subbc*gray_subbc);
real gray_subbjf = K(4)*(gray_subbjc*gray_subbjc);
real gray_subbjg = gray_subb*gray_subbaf + K(1);
real gray_subbjh = gray_subbdb*gray_subej + K(1);
real gray_subbji = gray_subbfd*gray_subej + K(1);
real gray_subbjj = gray_subbjg*gray_subbji;
real gray_subcaa = gray_subbjh*gray_subbjj;
real gray_subcab = gray_subej - K(1);
real gray_subcac = K(1.0)/(K(-48)*gray_subb*gray_subbb*gray_subbcg*gray_subbfc*gray_subbjb*K(1.0)/(gray_subbc*gray_subbc*gray_subbc*gray_subbc)*gray_subbjc*gray_subbjc*gray_subbjc*gray_subbjc + K(2)*gray_subb*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjd*gray_subbje*gray_subcab*gray_subdc + K(2)*gray_subb*gray_subbcg*gray_subbjd*gray_subbje*gray_subbji*gray_subdc*gray_subfc - gray_subb*gray_subbcg*gray_subbjf*gray_subbji*gray_subcab*gray_subcf*gray_subfc + K(2)*gray_subb*gray_subbfc*gray_subbjd*gray_subbje*gray_subbjh*gray_subdc*gray_subfc - gray_subb*gray_subbfc*gray_subbjf*gray_subbjh*gray_subcab*gray_subcf*gray_subfc - gray_subb*gray_subbjf*gray_subbjh*gray_subbji*gray_subcf - gray_subbb*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjf*gray_subbjg*gray_subcab*gray_subcf - gray_subbb*gray_subbcg*gray_subbjf*gray_subbjg*gray_subbji*gray_subcf*gray_subfc - gray_subbb*gray_subbfc*gray_subbjf*gray_subbjg*gray_subbjh*gray_subcf*gray_subfc + K(2)*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjd*gray_subbje*gray_subbjg*gray_subea + gray_subcaa*gray_subcab);
real gray_subcad = gray_subcac*(-gray_subb*gray_subbci*gray_subbfc*gray_subbjd*gray_subbje*gray_subdc*gray_subfi + gray_subb*gray_subbjf*gray_subbji*gray_subcf*gray_subeg*gray_subfi + gray_subbb*gray_subbci*gray_subbfc*gray_subbjf*gray_subbjg*gray_subcf*gray_subfi - gray_subbjj*gray_subfa*gray_subfi);
real gray_subcae = gray_subcac*(-gray_subb*gray_subbcg*gray_subbci*gray_subbjd*gray_subbje*gray_subdc*gray_subee + gray_subb*gray_subbjf*gray_subbjh*gray_subcf*gray_subee*gray_subeg + gray_subbb*gray_subbcg*gray_subbci*gray_subbjf*gray_subbjg*gray_subcf*gray_subee - gray_subbjg*gray_subbjh*gray_subee*gray_subfa);
real gray_subcaf = gray_subcac*(-gray_subbb*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjd*gray_subbje*z - gray_subbca*gray_subbjh*gray_subbji + gray_subbcg*gray_subbjf*gray_subbji*gray_subcf*gray_subdc*gray_subfc*z + gray_subbfc*gray_subbjf*gray_subbjh*gray_subcf*gray_subdc*gray_subfc*z);
real gray_subcag = gray_subcac*(-gray_subb*gray_subbjd*gray_subbje*gray_subdc*gray_subee*gray_subfc*gray_subfi + gray_subb*gray_subbjf*gray_subcab*gray_subcf*gray_subee*gray_subfc*gray_subfi + gray_subbb*gray_subbjf*gray_subbjg*gray_subcf*gray_subee*gray_subfc*gray_subfi - gray_subbdd*gray_subbjg*gray_subcab*gray_subee);
real gray_subcah = gray_subcac*(-gray_subbb*gray_subbci*gray_subbfc*gray_subbjd*gray_subbje*gray_subfi*z + gray_subbci*gray_subbfc*gray_subbjf*gray_subcab*gray_subcf*gray_subdc*gray_subfi*z - gray_subbgg*gray_subbji*gray_subcab*z + gray_subbjf*gray_subbji*gray_subcf*gray_subdc*gray_subeg*gray_subfi*z);
real gray_subcai = gray_subbjh*gray_subcab;
real gray_subcaj = gray_subcac*(-gray_subbb*gray_subbcg*gray_subbci*gray_subbjd*gray_subbje*gray_subee*z + gray_subbcg*gray_subbci*gray_subbjf*gray_subcab*gray_subcf*gray_subdc*gray_subee*z - gray_subbif*gray_subcai*z + gray_subbjf*gray_subbjh*gray_subcf*gray_subdc*gray_subee*gray_subeg*z);
real16 Xi = {K(0),
-gray_subdb - gray_subed - gray_subfh - gray_subgc,
-gray_subhe - gray_subib - gray_subih - gray_subja,
-gray_subbab - gray_subbai - gray_subbbf - gray_subbbh,
(K(1.0)/K(2.0))*gray_subdb + (K(1.0)/K(2.0))*gray_subed + (K(1.0)/K(2.0))*gray_subfh + (K(1.0)/K(2.0))*gray_subgc,
-gray_subbcd*gray_subbce - gray_subbcf*gray_subgb - gray_subbde*(gray_subbcg*gray_subbch + gray_subbcg*gray_subbda + gray_subbdb*gray_subdj - gray_subbdb*gray_subec + gray_subbdd*(gray_subbdc + gray_subbe*gray_subdi)) - gray_subbdi*gray_subbdj,
-gray_subbec + (K(1.0)/K(2.0))*gray_subbef - gray_subbeg + (K(1.0)/K(2.0))*gray_subbeh - gray_subbfb + (K(1.0)/K(2.0))*gray_subbfe - gray_subbfh + (K(1.0)/K(2.0))*gray_subbfi,
(K(1.0)/K(2.0))*gray_subbfj + (K(1.0)/K(2.0))*gray_subbgd + (K(1.0)/K(2.0))*gray_subbge + (K(1.0)/K(2.0))*gray_subbgf - gray_subbgj - gray_subbha - gray_subbhe - gray_subbhh,
(K(1.0)/K(2.0))*gray_subhe + (K(1.0)/K(2.0))*gray_subib + (K(1.0)/K(2.0))*gray_subih + (K(1.0)/K(2.0))*gray_subja,
(K(1.0)/K(2.0))*gray_subbec - gray_subbef + (K(1.0)/K(2.0))*gray_subbeg - gray_subbeh + (K(1.0)/K(2.0))*gray_subbfb - gray_subbfe + (K(1.0)/K(2.0))*gray_subbfh - gray_subbfi,
-gray_subbce*gray_subbhj - gray_subbcf*gray_subij - gray_subbde*gray_subbfg - gray_subbdj*(gray_subbdh*(gray_subbdc + gray_subgd*gray_subhi) + gray_subbei*gray_subbia + gray_subbfd*gray_subhf + gray_subbfd*gray_subhj - gray_subbfd*gray_subia),
(K(1.0)/K(2.0))*gray_subbib + (K(1.0)/K(2.0))*gray_subbic + (K(1.0)/K(2.0))*gray_subbid + (K(1.0)/K(2.0))*gray_subbie - gray_subbih - gray_subbii - gray_subbij - gray_subbja,
(K(1.0)/K(2.0))*gray_subbab + (K(1.0)/K(2.0))*gray_subbai + (K(1.0)/K(2.0))*gray_subbbf + (K(1.0)/K(2.0))*gray_subbbh,
-gray_subbfj - gray_subbgd - gray_subbge - gray_subbgf + (K(1.0)/K(2.0))*gray_subbgj + (K(1.0)/K(2.0))*gray_subbha + (K(1.0)/K(2.0))*gray_subbhe + (K(1.0)/K(2.0))*gray_subbhh,
-gray_subbib - gray_subbic - gray_subbid - gray_subbie + (K(1.0)/K(2.0))*gray_subbih + (K(1.0)/K(2.0))*gray_subbii + (K(1.0)/K(2.0))*gray_subbij + (K(1.0)/K(2.0))*gray_subbja,
-gray_subbaa*gray_subbcf - gray_subbce*(gray_subbaf*gray_subjb + gray_subbga*gray_subjf + gray_subbga*gray_subji + gray_subbgb*(-gray_subbad - gray_subjc)) - gray_subbde*gray_subbgi - gray_subbdj*gray_subbig};
real16 gUP = {gray_subcac*(K(2)*gray_subb*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjd*gray_subbje*gray_subdc - gray_subb*gray_subbcg*gray_subbjf*gray_subbji*gray_subcf*gray_subfc - gray_subb*gray_subbfc*gray_subbjf*gray_subbjh*gray_subcf*gray_subfc - gray_subbb*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjf*gray_subbjg*gray_subcf + gray_subcaa),
gray_subcad,
gray_subcae,
gray_subcaf,
gray_subcad,
gray_subcac*(K(2)*gray_subb*gray_subbfc*gray_subbjd*gray_subbje*gray_subdc*gray_subfc - gray_subb*gray_subbfc*gray_subbjf*gray_subcab*gray_subcf*gray_subfc - gray_subb*gray_subbjf*gray_subbji*gray_subcf - gray_subbb*gray_subbfc*gray_subbjf*gray_subbjg*gray_subcf*gray_subfc + gray_subbjj*gray_subcab),
gray_subcag,
gray_subcah,
gray_subcae,
gray_subcag,
gray_subcac*(K(2)*gray_subb*gray_subbcg*gray_subbjd*gray_subbje*gray_subdc*gray_subfc - gray_subb*gray_subbcg*gray_subbjf*gray_subcab*gray_subcf*gray_subfc - gray_subb*gray_subbjf*gray_subbjh*gray_subcf - gray_subbb*gray_subbcg*gray_subbjf*gray_subbjg*gray_subcf*gray_subfc + gray_subbjg*gray_subcai),
gray_subcaj,
gray_subcaf,
gray_subcah,
gray_subcaj,
gray_subcac*(-gray_subbb*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjf*gray_subcab*gray_subcf - gray_subbb*gray_subbcg*gray_subbjf*gray_subbji*gray_subcf*gray_subfc - gray_subbb*gray_subbfc*gray_subbjf*gray_subbjh*gray_subcf*gray_subfc + K(2)*gray_subbcg*gray_subbfc*gray_subbjb*gray_subbjd*gray_subbje*gray_subea + gray_subbji*gray_subcai)};
real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));
