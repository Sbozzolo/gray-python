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
real gray_subbc = GRAY_SQRT(gray_subbb);
real gray_subbd = gray_subc + gray_subbb*gray_subbb;
real gray_subbe = K(1.0)/gray_subbd;
real gray_subbf = K(2)*x;
real gray_subbg = K(1.0)*x;
real gray_subbh = K(1.0)/gray_subh;
real gray_subbi = gray_subbh*gray_subg;
real gray_subbj = gray_subbg*gray_subbi;
real gray_subca = gray_subbe*(gray_subbf + gray_subbj);
real gray_subcb = gray_subbc*gray_subca;
real gray_subcc = K(-1.0)*gray_suba + gray_subf + K(2)*gray_subh;
real gray_subcd = gray_subbe*gray_subcc;
real gray_subce = K(1.0)/gray_subbc;
real gray_subcf = K(0.5)*x;
real gray_subcg = gray_subbi*x;
real gray_subch = K(0.25)*gray_subcg;
real gray_subci = gray_subce*(gray_subcf + gray_subch);
real gray_subcj = gray_subcd*gray_subci;
real gray_subda = K(2.0)*x;
real gray_subdb = gray_subcg + gray_subda;
real gray_subdc = GRAY_SQRT_CUBE(gray_subbb);
real gray_subdd = K(1.0)/(gray_subbd*gray_subbd);
real gray_subde = gray_subcc*gray_subdd;
real gray_subdf = gray_subdc*gray_subde;
real gray_subdg = gray_subdb*gray_subdf;
real gray_subdh = a_spin*y + gray_subbc*x;
real gray_subdi = gray_subba + gray_subh + gray_subi;
real gray_subdj = K(1.0)/gray_subdi;
real gray_subea = gray_subcb*gray_subdj;
real gray_subeb = -gray_subbg - gray_subbi*gray_subcf;
real gray_subec = K(1.0)/(gray_subdi*gray_subdi);
real gray_subed = gray_subbc*gray_subcd;
real gray_subee = gray_subec*gray_subed;
real gray_subef = gray_subeb*gray_subee;
real gray_subeg = gray_subdh*gray_subdj;
real gray_subeh = gray_subbc + gray_subci*x;
real gray_subei = gray_subdj*gray_subed;
real gray_subej = gray_subcj*gray_subeg - gray_subdg*gray_subeg + gray_subdh*gray_subea + gray_subdh*gray_subef + gray_subeh*gray_subei;
real gray_subfa = -a_spin*x + gray_subbc*y;
real gray_subfb = -a_spin + gray_subci*y;
real gray_subfc = gray_subdj*gray_subfa;
real gray_subfd = gray_subcj*gray_subfc - gray_subdg*gray_subfc + gray_subea*gray_subfa + gray_subef*gray_subfa + gray_subei*gray_subfb;
real gray_subfe = gray_subca*z;
real gray_subff = gray_subdb*gray_subde;
real gray_subfg = gray_subbb*z;
real gray_subfh = gray_subff*gray_subfg;
real gray_subfi = gray_subfe - gray_subfh;
real gray_subfj = K(2)*y;
real gray_subga = K(1.0)*y;
real gray_subgb = gray_subbi*gray_subga;
real gray_subgc = gray_subbe*(gray_subfj + gray_subgb);
real gray_subgd = gray_subbc*gray_subgc;
real gray_subge = K(0.5)*y;
real gray_subgf = gray_subbi*y;
real gray_subgg = K(0.25)*gray_subgf;
real gray_subgh = gray_subce*(gray_subge + gray_subgg);
real gray_subgi = gray_subcd*gray_subgh;
real gray_subgj = K(2.0)*y;
real gray_subha = gray_subgf + gray_subgj;
real gray_subhb = gray_subdf*gray_subha;
real gray_subhc = gray_subdj*gray_subgd;
real gray_subhd = a_spin + gray_subgh*x;
real gray_subhe = -gray_subbi*gray_subge - gray_subga;
real gray_subhf = gray_subee*gray_subhe;
real gray_subhg = gray_subdh*gray_subhc + gray_subdh*gray_subhf + gray_subeg*gray_subgi - gray_subeg*gray_subhb + gray_subei*gray_subhd;
real gray_subhh = gray_subbc + gray_subgh*y;
real gray_subhi = gray_subei*gray_subhh + gray_subfa*gray_subhc + gray_subfa*gray_subhf + gray_subfc*gray_subgi - gray_subfc*gray_subhb;
real gray_subhj = gray_subgc*z;
real gray_subia = gray_subde*gray_subha;
real gray_subib = gray_subfg*gray_subia;
real gray_subic = gray_subhj - gray_subib;
real gray_subid = K(2)*z;
real gray_subie = K(0.5)*z;
real gray_subif = gray_subbh*(gray_suba*z + gray_subg*gray_subie);
real gray_subig = K(2)*gray_subif;
real gray_subih = gray_subbe*(gray_subid + gray_subig);
real gray_subii = gray_subbc*gray_subih;
real gray_subij = (K(1.0)/K(2.0))*gray_subif;
real gray_subja = gray_subie + gray_subij;
real gray_subjb = gray_subcd*gray_subce;
real gray_subjc = gray_subja*gray_subjb;
real gray_subjd = K(2.0)*z;
real gray_subje = gray_subde*(-gray_suba*gray_subid - gray_subbb*(gray_subig + gray_subjd));
real gray_subjf = gray_subbc*gray_subje;
real gray_subjg = gray_subcd*gray_subdj;
real gray_subjh = gray_subja*gray_subjg;
real gray_subji = gray_subdj*gray_subii;
real gray_subjj = -gray_subif - K(1.0)*z;
real gray_subbaa = gray_subee*gray_subjj;
real gray_subbab = gray_subbaa*gray_subdh + gray_subdh*gray_subji + gray_subeg*gray_subjc + gray_subeg*gray_subjf + gray_subjh*x;
real gray_subbac = gray_subbaa*gray_subfa + gray_subfa*gray_subji + gray_subfc*gray_subjc + gray_subfc*gray_subjf + gray_subjh*y;
real gray_subbad = gray_subih*z;
real gray_subbae = gray_subje*z;
real gray_subbaf = gray_subbad + gray_subbae + gray_subcd;
real gray_subbag = gray_subdh*gray_subdh;
real gray_subbah = gray_subbag*gray_subec;
real gray_subbai = -gray_subbj - gray_subda;
real gray_subbaj = K(1.0)/(gray_subdi*gray_subdi*gray_subdi);
real gray_subbba = gray_subbaj*gray_subed;
real gray_subbbb = gray_subbag*gray_subbba;
real gray_subbbc = K(2)*gray_subbc;
real gray_subbbd = gray_subdh*gray_subee;
real gray_subbbe = gray_subdh*gray_subec;
real gray_subbbf = gray_subbbe*gray_subfa;
real gray_subbbg = gray_subbba*gray_subdh*gray_subfa;
real gray_subbbh = gray_subee*gray_subfa;
real gray_subbbi = gray_subbai*gray_subbbg + gray_subbbd*gray_subfb + gray_subbbf*gray_subcb + gray_subbbf*gray_subcj - gray_subbbf*gray_subdg + gray_subbbh*gray_subeh;
real gray_subbbj = gray_subcd*z;
real gray_subbca = gray_subbbe*gray_subbbj;
real gray_subbcb = gray_subjg*z;
real gray_subbcc = gray_subbca*gray_subeb + gray_subbcb*gray_subeh + gray_subeg*gray_subfe - gray_subeg*gray_subfh;
real gray_subbcd = -gray_subgb - gray_subgj;
real gray_subbce = K(2)*a_spin;
real gray_subbcf = gray_subbbd*gray_subhh + gray_subbbf*gray_subgd + gray_subbbf*gray_subgi - gray_subbbf*gray_subhb + gray_subbbg*gray_subbcd + gray_subbbh*gray_subhd;
real gray_subbcg = gray_subbca*gray_subhe + gray_subbcb*gray_subhd + gray_subeg*gray_subhj - gray_subeg*gray_subib;
real gray_subbch = gray_subcd*gray_subja;
real gray_subbci = gray_subbbe*gray_subbch;
real gray_subbcj = -gray_subig - gray_subjd;
real gray_subbda = gray_subec*gray_subfa;
real gray_subbdb = gray_subbch*gray_subbda;
real gray_subbdc = gray_subbbf*gray_subii + gray_subbbf*gray_subjc + gray_subbbf*gray_subjf + gray_subbbg*gray_subbcj + gray_subbci*y + gray_subbdb*x;
real gray_subbdd = gray_subdh*gray_subjg;
real gray_subbde = gray_subdj*gray_subjc*z;
real gray_subbdf = gray_subbad*gray_subeg + gray_subbae*gray_subeg + gray_subbca*gray_subjj + gray_subbdd + gray_subbde*x;
real gray_subbdg = gray_subfa*gray_subfa;
real gray_subbdh = gray_subbdg*gray_subec;
real gray_subbdi = gray_subbba*gray_subbdg;
real gray_subbdj = gray_subbbj*gray_subbda;
real gray_subbea = gray_subbcb*gray_subfb + gray_subbdj*gray_subeb + gray_subfc*gray_subfe - gray_subfc*gray_subfh;
real gray_subbeb = gray_subbcb*gray_subhh + gray_subbdj*gray_subhe + gray_subfc*gray_subhj - gray_subfc*gray_subib;
real gray_subbec = gray_subfa*gray_subjg;
real gray_subbed = gray_subbad*gray_subfc + gray_subbae*gray_subfc + gray_subbde*y + gray_subbdj*gray_subjj + gray_subbec;
real gray_subbee = gray_subb*gray_subce;
real gray_subbef = gray_subb*gray_subcd*K(1.0)/gray_subdc;
real gray_subbeg = gray_subb*gray_subbc;
real gray_subbeh = K(1.0)/(gray_subdi*gray_subdi*gray_subdi*gray_subdi);
real gray_subbei = (K(1.0)/K(2.0))*gray_subb + (K(1.0)/K(2.0))*gray_subd + (K(1.0)/K(2.0))*gray_sube + gray_subj;
real gray_subbej = K(8)*(gray_subbei*gray_subbei*gray_subbei);
real gray_subbfa = K(1.0)/(gray_subbd*gray_subbd*gray_subbd);
real gray_subbfb = K(4)*(gray_subbei*gray_subbei);
real gray_subbfc = gray_subb*gray_subjb + K(1);
real gray_subbfd = gray_subbah*gray_subed + K(1);
real gray_subbfe = gray_subbdh*gray_subed + K(1);
real gray_subbff = gray_subbfc*gray_subbfe;
real gray_subbfg = gray_subbfd*gray_subbff;
real gray_subbfh = gray_subed - K(1);
real gray_subbfi = K(1.0)/(K(-48)*gray_subb*gray_subbag*gray_subbb*gray_subbdg*gray_subbeh*K(1.0)/(gray_subbd*gray_subbd*gray_subbd*gray_subbd)*gray_subbei*gray_subbei*gray_subbei*gray_subbei + K(2)*gray_subb*gray_subbag*gray_subbc*gray_subbdg*gray_subbeh*gray_subbej*gray_subbfa*gray_subbfh + K(2)*gray_subb*gray_subbag*gray_subbc*gray_subbej*gray_subbfa*gray_subbfe*gray_subec - gray_subb*gray_subbag*gray_subbfb*gray_subbfe*gray_subbfh*gray_subdd*gray_subec + K(2)*gray_subb*gray_subbc*gray_subbdg*gray_subbej*gray_subbfa*gray_subbfd*gray_subec - gray_subb*gray_subbdg*gray_subbfb*gray_subbfd*gray_subbfh*gray_subdd*gray_subec - gray_subb*gray_subbfb*gray_subbfd*gray_subbfe*gray_subdd - gray_subbag*gray_subbb*gray_subbdg*gray_subbeh*gray_subbfb*gray_subbfc*gray_subbfh*gray_subdd - gray_subbag*gray_subbb*gray_subbfb*gray_subbfc*gray_subbfe*gray_subdd*gray_subec + K(2)*gray_subbag*gray_subbdg*gray_subbeh*gray_subbej*gray_subbfa*gray_subbfc*gray_subdc - gray_subbb*gray_subbdg*gray_subbfb*gray_subbfc*gray_subbfd*gray_subdd*gray_subec + gray_subbfg*gray_subbfh);
real gray_subbfj = gray_subbfi*(-gray_subb*gray_subbaj*gray_subbc*gray_subbdg*gray_subbej*gray_subbfa*gray_subdh + gray_subb*gray_subbfb*gray_subbfe*gray_subdd*gray_subdh*gray_subdj + gray_subbaj*gray_subbb*gray_subbdg*gray_subbfb*gray_subbfc*gray_subdd*gray_subdh - gray_subbff*gray_subdh*gray_subei);
real gray_subbga = gray_subbfi*(-gray_subb*gray_subbag*gray_subbaj*gray_subbc*gray_subbej*gray_subbfa*gray_subfa + gray_subb*gray_subbfb*gray_subbfd*gray_subdd*gray_subdj*gray_subfa + gray_subbag*gray_subbaj*gray_subbb*gray_subbfb*gray_subbfc*gray_subdd*gray_subfa - gray_subbfc*gray_subbfd*gray_subei*gray_subfa);
real gray_subbgb = gray_subbfi*(-gray_subbag*gray_subbb*gray_subbdg*gray_subbeh*gray_subbej*gray_subbfa*z + gray_subbag*gray_subbc*gray_subbfb*gray_subbfe*gray_subdd*gray_subec*z - gray_subbbj*gray_subbfd*gray_subbfe + gray_subbc*gray_subbdg*gray_subbfb*gray_subbfd*gray_subdd*gray_subec*z);
real gray_subbgc = gray_subbfi*(-gray_subb*gray_subbc*gray_subbej*gray_subbfa*gray_subdh*gray_subec*gray_subfa + gray_subb*gray_subbfb*gray_subbfh*gray_subdd*gray_subdh*gray_subec*gray_subfa + gray_subbb*gray_subbfb*gray_subbfc*gray_subdd*gray_subdh*gray_subec*gray_subfa - gray_subbbd*gray_subbfc*gray_subbfh*gray_subfa);
real gray_subbgd = gray_subbfi*(-gray_subbaj*gray_subbb*gray_subbdg*gray_subbej*gray_subbfa*gray_subdh*z + gray_subbaj*gray_subbc*gray_subbdg*gray_subbfb*gray_subbfh*gray_subdd*gray_subdh*z + gray_subbc*gray_subbfb*gray_subbfe*gray_subdd*gray_subdh*gray_subdj*z - gray_subbdd*gray_subbfe*gray_subbfh*z);
real gray_subbge = gray_subbfd*gray_subbfh;
real gray_subbgf = gray_subbfi*(-gray_subbag*gray_subbaj*gray_subbb*gray_subbej*gray_subbfa*gray_subfa*z + gray_subbag*gray_subbaj*gray_subbc*gray_subbfb*gray_subbfh*gray_subdd*gray_subfa*z + gray_subbc*gray_subbfb*gray_subbfd*gray_subdd*gray_subdj*gray_subfa*z - gray_subbec*gray_subbge*z);
real16 Phi_t = {K(0),
K(0),
K(0),
K(0),
gray_subcb + gray_subcj - gray_subdg,
gray_subej,
gray_subfd,
gray_subfi,
gray_subgd + gray_subgi - gray_subhb,
gray_subhg,
gray_subhi,
gray_subic,
gray_subii + gray_subjc + gray_subjf,
gray_subbab,
gray_subbac,
gray_subbaf};
real16 Phi_x = {K(0),
K(0),
K(0),
K(0),
gray_subej,
gray_subbah*gray_subcb + gray_subbah*gray_subcj - gray_subbah*gray_subdg + gray_subbai*gray_subbbb + gray_subbbd*(gray_subbbc + gray_subbf*gray_subci),
gray_subbbi,
gray_subbcc,
gray_subhg,
gray_subbah*gray_subgd + gray_subbah*gray_subgi - gray_subbah*gray_subhb + gray_subbbb*gray_subbcd + gray_subbbd*(gray_subbce + gray_subbf*gray_subgh),
gray_subbcf,
gray_subbcg,
gray_subbab,
gray_subbah*gray_subii + gray_subbah*gray_subjc + gray_subbah*gray_subjf + gray_subbbb*gray_subbcj + gray_subbci*gray_subbf,
gray_subbdc,
gray_subbdf};
real16 Phi_y = {K(0),
K(0),
K(0),
K(0),
gray_subfd,
gray_subbbi,
gray_subbai*gray_subbdi + gray_subbbh*(-gray_subbce + gray_subci*gray_subfj) + gray_subbdh*gray_subcb + gray_subbdh*gray_subcj - gray_subbdh*gray_subdg,
gray_subbea,
gray_subhi,
gray_subbcf,
gray_subbbh*(gray_subbbc + gray_subfj*gray_subgh) + gray_subbcd*gray_subbdi + gray_subbdh*gray_subgd + gray_subbdh*gray_subgi - gray_subbdh*gray_subhb,
gray_subbeb,
gray_subbac,
gray_subbdc,
gray_subbcj*gray_subbdi + gray_subbdb*gray_subfj + gray_subbdh*gray_subii + gray_subbdh*gray_subjc + gray_subbdh*gray_subjf,
gray_subbed};
real16 Phi_z = {K(0),
K(0),
K(0),
K(0),
gray_subfi,
gray_subbcc,
gray_subbea,
gray_subbee*gray_subca + gray_subbef*(-gray_subcf - gray_subch) - gray_subbeg*gray_subff,
gray_subic,
gray_subbcg,
gray_subbeb,
gray_subbee*gray_subgc + gray_subbef*(-gray_subge - gray_subgg) - gray_subbeg*gray_subia,
gray_subbaf,
gray_subbdf,
gray_subbed,
gray_subbee*gray_subih + gray_subbee*gray_subje + gray_subbef*(-gray_subie - gray_subij) + gray_subid*gray_subjb};
real16 gUP = {gray_subbfi*(K(2)*gray_subb*gray_subbag*gray_subbc*gray_subbdg*gray_subbeh*gray_subbej*gray_subbfa - gray_subb*gray_subbag*gray_subbfb*gray_subbfe*gray_subdd*gray_subec - gray_subb*gray_subbdg*gray_subbfb*gray_subbfd*gray_subdd*gray_subec - gray_subbag*gray_subbb*gray_subbdg*gray_subbeh*gray_subbfb*gray_subbfc*gray_subdd + gray_subbfg),
gray_subbfj,
gray_subbga,
gray_subbgb,
gray_subbfj,
gray_subbfi*(K(2)*gray_subb*gray_subbc*gray_subbdg*gray_subbej*gray_subbfa*gray_subec - gray_subb*gray_subbdg*gray_subbfb*gray_subbfh*gray_subdd*gray_subec - gray_subb*gray_subbfb*gray_subbfe*gray_subdd - gray_subbb*gray_subbdg*gray_subbfb*gray_subbfc*gray_subdd*gray_subec + gray_subbff*gray_subbfh),
gray_subbgc,
gray_subbgd,
gray_subbga,
gray_subbgc,
gray_subbfi*(K(2)*gray_subb*gray_subbag*gray_subbc*gray_subbej*gray_subbfa*gray_subec - gray_subb*gray_subbag*gray_subbfb*gray_subbfh*gray_subdd*gray_subec - gray_subb*gray_subbfb*gray_subbfd*gray_subdd - gray_subbag*gray_subbb*gray_subbfb*gray_subbfc*gray_subdd*gray_subec + gray_subbfc*gray_subbge),
gray_subbgf,
gray_subbgb,
gray_subbgd,
gray_subbgf,
gray_subbfi*(-gray_subbag*gray_subbb*gray_subbdg*gray_subbeh*gray_subbfb*gray_subbfh*gray_subdd - gray_subbag*gray_subbb*gray_subbfb*gray_subbfe*gray_subdd*gray_subec + K(2)*gray_subbag*gray_subbdg*gray_subbeh*gray_subbej*gray_subbfa*gray_subdc - gray_subbb*gray_subbdg*gray_subbfb*gray_subbfd*gray_subdd*gray_subec + gray_subbfe*gray_subbge)};
real16 Upsilon = -uUPt * Phi_t - uUPx * Phi_x - uUPy * Phi_y - uUPz * Phi_z;
real16 UpsilonT = Upsilon.s048c159d26ae37bf;
real16 Xi = UpsilonT - Upsilon/2;
real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));
