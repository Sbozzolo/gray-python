real gray_suba = z*z;
real gray_subb = a_spin*a_spin;
real gray_subc = gray_suba*gray_subb;
real gray_subd = x*x;
real gray_sube = y*y;
real gray_subf = gray_suba + gray_subd + gray_sube;
real gray_subg = -gray_subb + gray_subf;
real gray_subh = GRAY_SQRT(gray_subc + K(0.25)*(gray_subg*gray_subg));
real gray_subi = K(0.5)*gray_subb;
real gray_subj = gray_subh - gray_subi;
real gray_subba = K(0.5)*gray_suba + K(0.5)*gray_subd + K(0.5)*gray_sube;
real gray_subbb = gray_subba + gray_subj;
real gray_subbc = GRAY_SQRT(gray_subbb);
real gray_subbd = gray_subba + gray_subh + gray_subi;
real gray_subbe = K(1.0)/(gray_subbd*gray_subbd*gray_subbd);
real gray_subbf = (K(1.0)/K(2.0))*gray_suba + (K(1.0)/K(2.0))*gray_subd + (K(1.0)/K(2.0))*gray_sube + gray_subj;
real gray_subbg = K(8)*(gray_subbf*gray_subbf*gray_subbf);
real gray_subbh = a_spin*y + gray_subbc*x;
real gray_subbi = -a_spin*x + gray_subbc*y;
real gray_subbj = gray_subbi*gray_subbi;
real gray_subca = gray_subc + gray_subbb*gray_subbb;
real gray_subcb = K(1.0)/(gray_subca*gray_subca*gray_subca);
real gray_subcc = K(4)*(gray_subbf*gray_subbf);
real gray_subcd = K(1.0)/(gray_subca*gray_subca);
real gray_subce = K(1.0)/gray_subbc;
real gray_subcf = K(-1.0)*gray_subb + gray_subf + K(2)*gray_subh;
real gray_subcg = K(1.0)/gray_subca;
real gray_subch = gray_subcf*gray_subcg;
real gray_subci = gray_subce*gray_subch;
real gray_subcj = gray_suba*gray_subci + K(1);
real gray_subda = K(1.0)/gray_subbd;
real gray_subdb = K(1.0)/(gray_subbd*gray_subbd);
real gray_subdc = gray_subbc*gray_subch;
real gray_subdd = gray_subdb*gray_subdc;
real gray_subde = gray_subbj*gray_subdd + K(1);
real gray_subdf = gray_subcj*gray_subde;
real gray_subdg = gray_subda*gray_subdc;
real gray_subdh = -gray_suba*gray_subbc*gray_subbe*gray_subbg*gray_subbh*gray_subbj*gray_subcb + gray_suba*gray_subbh*gray_subcc*gray_subcd*gray_subda*gray_subde + gray_subbb*gray_subbe*gray_subbh*gray_subbj*gray_subcc*gray_subcd*gray_subcj - gray_subbh*gray_subdf*gray_subdg;
real gray_subdi = K(2)*x;
real gray_subdj = K(1.0)*x;
real gray_subea = K(1.0)/gray_subh;
real gray_subeb = gray_subea*gray_subg;
real gray_subec = gray_subdj*gray_subeb;
real gray_subed = gray_subdi + gray_subec;
real gray_subee = gray_subbc*gray_subcg;
real gray_subef = gray_subed*gray_subee;
real gray_subeg = K(0.5)*x;
real gray_subeh = gray_subeb*x;
real gray_subei = K(0.25)*gray_subeh;
real gray_subej = gray_subce*(gray_subeg + gray_subei);
real gray_subfa = gray_subch*gray_subej;
real gray_subfb = K(2.0)*x;
real gray_subfc = gray_subeh + gray_subfb;
real gray_subfd = GRAY_SQRT_CUBE(gray_subbb);
real gray_subfe = gray_subcd*gray_subcf;
real gray_subff = gray_subfd*gray_subfe;
real gray_subfg = gray_subfc*gray_subff;
real gray_subfh = K(1.0)/(gray_subbd*gray_subbd*gray_subbd*gray_subbd);
real gray_subfi = gray_subbh*gray_subbh;
real gray_subfj = gray_subdc - K(1);
real gray_subga = gray_subdd*gray_subfi + K(1);
real gray_subgb = gray_subdf*gray_subga;
real gray_subgc = K(1.0)/(K(-48)*gray_suba*gray_subbb*gray_subbj*gray_subfh*gray_subfi*gray_subbf*gray_subbf*gray_subbf*gray_subbf*K(1.0)/(gray_subca*gray_subca*gray_subca*gray_subca) + K(2)*gray_suba*gray_subbc*gray_subbg*gray_subbj*gray_subcb*gray_subdb*gray_subga + K(2)*gray_suba*gray_subbc*gray_subbg*gray_subbj*gray_subcb*gray_subfh*gray_subfi*gray_subfj + K(2)*gray_suba*gray_subbc*gray_subbg*gray_subcb*gray_subdb*gray_subde*gray_subfi - gray_suba*gray_subbj*gray_subcc*gray_subcd*gray_subdb*gray_subfj*gray_subga - gray_suba*gray_subcc*gray_subcd*gray_subdb*gray_subde*gray_subfi*gray_subfj - gray_suba*gray_subcc*gray_subcd*gray_subde*gray_subga - gray_subbb*gray_subbj*gray_subcc*gray_subcd*gray_subcj*gray_subdb*gray_subga - gray_subbb*gray_subbj*gray_subcc*gray_subcd*gray_subcj*gray_subfh*gray_subfi*gray_subfj - gray_subbb*gray_subcc*gray_subcd*gray_subcj*gray_subdb*gray_subde*gray_subfi + K(2)*gray_subbg*gray_subbj*gray_subcb*gray_subcj*gray_subfd*gray_subfh*gray_subfi + gray_subfj*gray_subgb);
real gray_subgd = K(0.5)*gray_subgc;
real gray_subge = gray_subgd*(gray_subef + gray_subfa - gray_subfg);
real gray_subgf = gray_subdh*gray_subge;
real gray_subgg = K(2)*y;
real gray_subgh = K(1.0)*y;
real gray_subgi = gray_subeb*gray_subgh;
real gray_subgj = gray_subgg + gray_subgi;
real gray_subha = gray_subee*gray_subgj;
real gray_subhb = K(0.5)*y;
real gray_subhc = gray_subeb*y;
real gray_subhd = K(0.25)*gray_subhc;
real gray_subhe = gray_subce*(gray_subhb + gray_subhd);
real gray_subhf = gray_subch*gray_subhe;
real gray_subhg = K(2.0)*y;
real gray_subhh = gray_subhc + gray_subhg;
real gray_subhi = gray_subff*gray_subhh;
real gray_subhj = gray_subha + gray_subhf - gray_subhi;
real gray_subia = -gray_suba*gray_subbc*gray_subbe*gray_subbg*gray_subbi*gray_subcb*gray_subfi + gray_suba*gray_subbi*gray_subcc*gray_subcd*gray_subda*gray_subga + gray_subbb*gray_subbe*gray_subbi*gray_subcc*gray_subcd*gray_subcj*gray_subfi - gray_subbi*gray_subcj*gray_subdg*gray_subga;
real gray_subib = gray_subgd*gray_subia;
real gray_subic = gray_subhj*gray_subib;
real gray_subid = K(2)*z;
real gray_subie = K(0.5)*z;
real gray_subif = gray_subea*(gray_subb*z + gray_subg*gray_subie);
real gray_subig = K(2)*gray_subif;
real gray_subih = gray_subid + gray_subig;
real gray_subii = gray_subee*gray_subih;
real gray_subij = (K(1.0)/K(2.0))*gray_subif;
real gray_subja = gray_subie + gray_subij;
real gray_subjb = gray_subci*gray_subja;
real gray_subjc = K(2.0)*z;
real gray_subjd = gray_subfe*(-gray_subb*gray_subid - gray_subbb*(gray_subig + gray_subjc));
real gray_subje = gray_subbc*gray_subjd;
real gray_subjf = gray_subii + gray_subjb + gray_subje;
real gray_subjg = gray_subch*z;
real gray_subjh = -gray_subbb*gray_subbg*gray_subbj*gray_subcb*gray_subfh*gray_subfi*z + gray_subbc*gray_subbj*gray_subcc*gray_subcd*gray_subdb*gray_subga*z + gray_subbc*gray_subcc*gray_subcd*gray_subdb*gray_subde*gray_subfi*z - gray_subde*gray_subga*gray_subjg;
real gray_subji = gray_subgd*gray_subjh;
real gray_subjj = gray_subjf*gray_subji;
real gray_subbaa = gray_subcg*z;
real gray_subbab = gray_subbaa*gray_subed;
real gray_subbac = gray_subfc*gray_subfe;
real gray_subbad = gray_subbb*z;
real gray_subbae = gray_subbac*gray_subbad;
real gray_subbaf = gray_subbab - gray_subbae;
real gray_subbag = gray_subbaf*gray_subji;
real gray_subbah = K(2)*gray_suba*gray_subbc*gray_subbg*gray_subbj*gray_subcb*gray_subfh*gray_subfi - gray_suba*gray_subbj*gray_subcc*gray_subcd*gray_subdb*gray_subga - gray_suba*gray_subcc*gray_subcd*gray_subdb*gray_subde*gray_subfi - gray_subbb*gray_subbj*gray_subcc*gray_subcd*gray_subcj*gray_subfh*gray_subfi + gray_subgb;
real gray_subbai = gray_subbh*gray_subda;
real gray_subbaj = a_spin + gray_subhe*x;
real gray_subbba = -gray_subeb*gray_subhb - gray_subgh;
real gray_subbbb = gray_subbh*gray_subdd;
real gray_subbbc = gray_subbai*gray_subha + gray_subbai*gray_subhf - gray_subbai*gray_subhi + gray_subbaj*gray_subdg + gray_subbba*gray_subbbb;
real gray_subbbd = gray_subbbc*gray_subib;
real gray_subbbe = gray_subbi*gray_subda;
real gray_subbbf = -a_spin + gray_subej*y;
real gray_subbbg = -gray_subdj - gray_subeb*gray_subeg;
real gray_subbbh = gray_subbi*gray_subdd;
real gray_subbbi = gray_subbbe*gray_subef + gray_subbbe*gray_subfa - gray_subbbe*gray_subfg + gray_subbbf*gray_subdg + gray_subbbg*gray_subbbh;
real gray_subbbj = gray_subbbi*gray_subib;
real gray_subbca = gray_subch*gray_subda;
real gray_subbcb = gray_subbca*gray_subja;
real gray_subbcc = -gray_subif - K(1.0)*z;
real gray_subbcd = gray_subbai*gray_subii + gray_subbai*gray_subjb + gray_subbai*gray_subje + gray_subbbb*gray_subbcc + gray_subbcb*x;
real gray_subbce = gray_subbcd*gray_subji;
real gray_subbcf = gray_subbag + gray_subbah*gray_subge - gray_subbbd + gray_subbbj - gray_subbce;
real gray_subbcg = gray_subbaa*gray_subgj;
real gray_subbch = gray_subfe*gray_subhh;
real gray_subbci = gray_subbad*gray_subbch;
real gray_subbcj = gray_subbcg - gray_subbci;
real gray_subbda = gray_subbcj*gray_subji;
real gray_subbdb = gray_subbah*gray_subgd;
real gray_subbdc = gray_subdh*gray_subgd;
real gray_subbdd = gray_subbbc*gray_subbdc;
real gray_subbde = gray_subbbi*gray_subbdc;
real gray_subbdf = gray_subbbe*gray_subii + gray_subbbe*gray_subjb + gray_subbbe*gray_subje + gray_subbbh*gray_subbcc + gray_subbcb*y;
real gray_subbdg = gray_subbdf*gray_subji;
real gray_subbdh = gray_subbda + gray_subbdb*gray_subhj + gray_subbdd - gray_subbde - gray_subbdg;
real gray_subbdi = gray_subbcj*gray_subib;
real gray_subbdj = gray_subbaf*gray_subbdc;
real gray_subbea = gray_subbcd*gray_subbdc;
real gray_subbeb = gray_subbdf*gray_subib;
real gray_subbec = gray_subbdb*gray_subjf - gray_subbdi - gray_subbdj + gray_subbea + gray_subbeb;
real gray_subbed = gray_subbh*gray_subdb;
real gray_subbee = gray_subbed*gray_subjg;
real gray_subbef = gray_subbc + gray_subej*x;
real gray_subbeg = gray_subbca*z;
real gray_subbeh = gray_subbab*gray_subbai - gray_subbae*gray_subbai + gray_subbbg*gray_subbee + gray_subbef*gray_subbeg;
real gray_subbei = K(1.0)*gray_subgc;
real gray_subbej = gray_subbei*gray_subjh;
real gray_subbfa = gray_subdb*gray_subfi;
real gray_subbfb = -gray_subgi - gray_subhg;
real gray_subbfc = gray_subbe*gray_subdc;
real gray_subbfd = gray_subbfc*gray_subfi;
real gray_subbfe = K(2)*a_spin;
real gray_subbff = gray_subbbb*(gray_subbfe + gray_subdi*gray_subhe) + gray_subbfa*gray_subha + gray_subbfa*gray_subhf - gray_subbfa*gray_subhi + gray_subbfb*gray_subbfd;
real gray_subbfg = gray_subch*gray_subja;
real gray_subbfh = gray_subbed*gray_subbfg;
real gray_subbfi = -gray_subig - gray_subjc;
real gray_subbfj = gray_subbfa*gray_subii + gray_subbfa*gray_subjb + gray_subbfa*gray_subje + gray_subbfd*gray_subbfi + gray_subbfh*gray_subdi;
real gray_subbga = -gray_subec - gray_subfb;
real gray_subbgb = K(2)*gray_subbc;
real gray_subbgc = gray_subbbb*(gray_subbgb + gray_subdi*gray_subej) + gray_subbfa*gray_subef + gray_subbfa*gray_subfa - gray_subbfa*gray_subfg + gray_subbfd*gray_subbga;
real gray_subbgd = gray_subbai*gray_subef + gray_subbai*gray_subfa - gray_subbai*gray_subfg + gray_subbbb*gray_subbbg + gray_subbef*gray_subdg;
real gray_subbge = gray_subbah*gray_subbei;
real gray_subbgf = gray_subbed*gray_subbi;
real gray_subbgg = gray_subbfc*gray_subbh*gray_subbi;
real gray_subbgh = gray_subbbb*gray_subbbf + gray_subbbh*gray_subbef + gray_subbga*gray_subbgg + gray_subbgf*gray_subef + gray_subbgf*gray_subfa - gray_subbgf*gray_subfg;
real gray_subbgi = gray_subbei*gray_subia;
real gray_subbgj = gray_subbai*gray_subbcg - gray_subbai*gray_subbci + gray_subbaj*gray_subbeg + gray_subbba*gray_subbee;
real gray_subbha = gray_subbi*gray_subdb;
real gray_subbhb = gray_subbha*gray_subjg;
real gray_subbhc = gray_subbab*gray_subbbe - gray_subbae*gray_subbbe + gray_subbbf*gray_subbeg + gray_subbbg*gray_subbhb;
real gray_subbhd = gray_subbj*gray_subdb;
real gray_subbhe = gray_subbfc*gray_subbj;
real gray_subbhf = gray_subbbh*(-gray_subbfe + gray_subej*gray_subgg) + gray_subbga*gray_subbhe + gray_subbhd*gray_subef + gray_subbhd*gray_subfa - gray_subbhd*gray_subfg;
real gray_subbhg = gray_subbfg*gray_subbha;
real gray_subbhh = gray_subbfh*y + gray_subbfi*gray_subbgg + gray_subbgf*gray_subii + gray_subbgf*gray_subjb + gray_subbgf*gray_subje + gray_subbhg*x;
real gray_subbhi = gray_subbbc*gray_subbdb + gray_subbbi*gray_subbdb + gray_subbdc*gray_subbff + gray_subbgj*gray_subji + gray_subbhc*gray_subji + gray_subbhf*gray_subib - gray_subbhh*gray_subji;
real gray_subbhj = gray_suba*gray_subce;
real gray_subbia = gray_subbhj*gray_subcg;
real gray_subbib = gray_suba*gray_subch*K(1.0)/gray_subfd;
real gray_subbic = gray_suba*gray_subbc;
real gray_subbid = -gray_subbac*gray_subbic + gray_subbia*gray_subed + gray_subbib*(-gray_subeg - gray_subei);
real gray_subbie = gray_subbaf*gray_subbdb + gray_subbcd*gray_subbdb + gray_subbdc*gray_subbfj - gray_subbgj*gray_subib + gray_subbhc*gray_subib + gray_subbhh*gray_subib + gray_subbid*gray_subji;
real gray_subbif = gray_subbc + gray_subhe*y;
real gray_subbig = gray_subbba*gray_subbhb + gray_subbbe*gray_subbcg - gray_subbbe*gray_subbci + gray_subbeg*gray_subbif;
real gray_subbih = gray_subbfi*gray_subbhe + gray_subbhd*gray_subii + gray_subbhd*gray_subjb + gray_subbhd*gray_subje + gray_subbhg*gray_subgg;
real gray_subbii = gray_subbbh*(gray_subbgb + gray_subgg*gray_subhe) + gray_subbfb*gray_subbhe + gray_subbhd*gray_subha + gray_subbhd*gray_subhf - gray_subbhd*gray_subhi;
real gray_subbij = gray_subbba*gray_subbbh + gray_subbbe*gray_subha + gray_subbbe*gray_subhf - gray_subbbe*gray_subhi + gray_subbif*gray_subdg;
real gray_subbja = gray_subbaj*gray_subbbh + gray_subbbb*gray_subbif + gray_subbfb*gray_subbgg + gray_subbgf*gray_subha + gray_subbgf*gray_subhf - gray_subbgf*gray_subhi;
real gray_subbjb = gray_subbei*gray_subdh;
real gray_subbjc = -gray_subbch*gray_subbic + gray_subbia*gray_subgj + gray_subbib*(-gray_subhb - gray_subhd);
real gray_subbjd = gray_subbcj*gray_subbdb + gray_subbdb*gray_subbdf + gray_subbdc*gray_subbgj - gray_subbdc*gray_subbhc + gray_subbdc*gray_subbhh + gray_subbih*gray_subib + gray_subbjc*gray_subji;
real gray_subbje = gray_subbhj*gray_subjd + gray_subbia*gray_subih + gray_subbib*(-gray_subie - gray_subij) + gray_subci*gray_subid;
real gray_subbjf = gray_subbaa*gray_subih;
real gray_subbjg = gray_subjd*z;
real gray_subbjh = gray_subbjf + gray_subbjg + gray_subch;
real gray_subbji = gray_subbca*gray_subbh;
real gray_subbjj = gray_subda*gray_subjb*z;
real gray_subcaa = gray_subbai*gray_subbjf + gray_subbai*gray_subbjg + gray_subbcc*gray_subbee + gray_subbji + gray_subbjj*x;
real gray_subcab = gray_subbca*gray_subbi;
real gray_subcac = gray_subbbe*gray_subbjf + gray_subbbe*gray_subbjg + gray_subbcc*gray_subbhb + gray_subbjj*y + gray_subcab;
real gray_subcad = -gray_suba*gray_subbc*gray_subbg*gray_subbh*gray_subbi*gray_subcb*gray_subdb + gray_suba*gray_subbh*gray_subbi*gray_subcc*gray_subcd*gray_subdb*gray_subfj + gray_subbb*gray_subbh*gray_subbi*gray_subcc*gray_subcd*gray_subcj*gray_subdb - gray_subbbh*gray_subbh*gray_subcj*gray_subfj;
real gray_subcae = gray_subcad*gray_subgd;
real gray_subcaf = K(2)*gray_suba*gray_subbc*gray_subbg*gray_subbj*gray_subcb*gray_subdb - gray_suba*gray_subbj*gray_subcc*gray_subcd*gray_subdb*gray_subfj - gray_suba*gray_subcc*gray_subcd*gray_subde - gray_subbb*gray_subbj*gray_subcc*gray_subcd*gray_subcj*gray_subdb + gray_subdf*gray_subfj;
real gray_subcag = -gray_subbb*gray_subbe*gray_subbg*gray_subbh*gray_subbj*gray_subcb*z + gray_subbc*gray_subbe*gray_subbh*gray_subbj*gray_subcc*gray_subcd*gray_subfj*z + gray_subbc*gray_subbh*gray_subcc*gray_subcd*gray_subda*gray_subde*z - gray_subbji*gray_subde*gray_subfj*z;
real gray_subcah = gray_subcag*gray_subgd;
real gray_subcai = gray_subbaf*gray_subcah;
real gray_subcaj = gray_subbbc*gray_subcae;
real gray_subcba = gray_subbbi*gray_subcae;
real gray_subcbb = gray_subbcd*gray_subcah;
real gray_subcbc = gray_subcai - gray_subcaj + gray_subcba - gray_subcbb + gray_subgf;
real gray_subcbd = gray_subcaf*gray_subgd;
real gray_subcbe = gray_subbbc*gray_subcbd - gray_subbbi*gray_subcbd + gray_subbcj*gray_subcah + gray_subbdc*gray_subhj - gray_subbdf*gray_subcah;
real gray_subcbf = -gray_subbaf*gray_subcbd + gray_subbcd*gray_subcbd - gray_subbcj*gray_subcae + gray_subbdc*gray_subjf + gray_subbdf*gray_subcae;
real gray_subcbg = gray_subbei*gray_subcag;
real gray_subcbh = gray_subbff*gray_subcae;
real gray_subcbi = gray_subbfj*gray_subcah;
real gray_subcbj = gray_subbei*gray_subcad;
real gray_subcca = gray_subbgj*gray_subcah;
real gray_subccb = gray_subbhc*gray_subcah;
real gray_subccc = gray_subbhf*gray_subcae;
real gray_subccd = gray_subbhh*gray_subcah;
real gray_subcce = gray_subbdd + gray_subbde + gray_subbff*gray_subcbd + gray_subcca + gray_subccb + gray_subccc - gray_subccd;
real gray_subccf = gray_subbid*gray_subcah;
real gray_subccg = gray_subbgj*gray_subcae;
real gray_subcch = gray_subbhc*gray_subcae;
real gray_subcci = gray_subbhh*gray_subcae;
real gray_subccj = gray_subbdj + gray_subbea + gray_subbfj*gray_subcbd + gray_subccf - gray_subccg + gray_subcch + gray_subcci;
real gray_subcda = gray_subbei*gray_subcaf;
real gray_subcdb = gray_subbcj*gray_subbdc + gray_subbdc*gray_subbdf + gray_subbgj*gray_subcbd - gray_subbhc*gray_subcbd + gray_subbhh*gray_subcbd + gray_subbih*gray_subcae + gray_subbjc*gray_subcah;
real gray_subcdc = gray_subfj*gray_subga;
real gray_subcdd = K(2)*gray_suba*gray_subbc*gray_subbg*gray_subcb*gray_subdb*gray_subfi - gray_suba*gray_subcc*gray_subcd*gray_subdb*gray_subfi*gray_subfj - gray_suba*gray_subcc*gray_subcd*gray_subga - gray_subbb*gray_subcc*gray_subcd*gray_subcj*gray_subdb*gray_subfi + gray_subcdc*gray_subcj;
real gray_subcde = gray_subcdd*gray_subgd;
real gray_subcdf = -gray_subbb*gray_subbe*gray_subbg*gray_subbi*gray_subcb*gray_subfi*z + gray_subbc*gray_subbe*gray_subbi*gray_subcc*gray_subcd*gray_subfi*gray_subfj*z + gray_subbc*gray_subbi*gray_subcc*gray_subcd*gray_subda*gray_subga*z - gray_subcab*gray_subcdc*z;
real gray_subcdg = gray_subcdf*gray_subgd;
real gray_subcdh = gray_subbaf*gray_subcdg - gray_subbbc*gray_subcde + gray_subbbi*gray_subcde - gray_subbcd*gray_subcdg + gray_subge*gray_subia;
real gray_subcdi = gray_subbcj*gray_subcdg;
real gray_subcdj = gray_subbdf*gray_subcdg;
real gray_subcea = gray_subcaj - gray_subcba + gray_subcdi - gray_subcdj + gray_subic;
real gray_subceb = -gray_subbaf*gray_subcae + gray_subbcd*gray_subcae - gray_subbcj*gray_subcde + gray_subbdf*gray_subcde + gray_subib*gray_subjf;
real gray_subcec = gray_subbei*gray_subcdf;
real gray_subced = gray_subbei*gray_subcdd;
real gray_subcee = gray_subbgj*gray_subcdg;
real gray_subcef = gray_subbhc*gray_subcdg;
real gray_subceg = gray_subbhh*gray_subcdg;
real gray_subceh = gray_subbbd + gray_subbbj + gray_subbhf*gray_subcde + gray_subcbh + gray_subcee + gray_subcef - gray_subceg;
real gray_subcei = gray_subbaf*gray_subib + gray_subbcd*gray_subib + gray_subbfj*gray_subcae - gray_subbgj*gray_subcde + gray_subbhc*gray_subcde + gray_subbhh*gray_subcde + gray_subbid*gray_subcdg;
real gray_subcej = gray_subbih*gray_subcdg;
real gray_subcfa = gray_subbjc*gray_subcdg;
real gray_subcfb = gray_subbdi + gray_subbeb + gray_subbih*gray_subcde + gray_subccg - gray_subcch + gray_subcci + gray_subcfa;
real gray_subcfc = -gray_subbb*gray_subbj*gray_subcc*gray_subcd*gray_subdb*gray_subga - gray_subbb*gray_subbj*gray_subcc*gray_subcd*gray_subfh*gray_subfi*gray_subfj - gray_subbb*gray_subcc*gray_subcd*gray_subdb*gray_subde*gray_subfi + K(2)*gray_subbg*gray_subbj*gray_subcb*gray_subfd*gray_subfh*gray_subfi + gray_subcdc*gray_subde;
real gray_subcfd = gray_subcfc*gray_subgd;
real gray_subcfe = gray_subbaf*gray_subcfd - gray_subbbc*gray_subcdg + gray_subbbi*gray_subcdg - gray_subbcd*gray_subcfd + gray_subge*gray_subjh;
real gray_subcff = gray_subbbc*gray_subcah - gray_subbbi*gray_subcah + gray_subbcj*gray_subcfd - gray_subbdf*gray_subcfd + gray_subhj*gray_subji;
real gray_subcfg = -gray_subcai + gray_subcbb - gray_subcdi + gray_subcdj + gray_subjj;
real gray_subcfh = gray_subbei*gray_subcfc;
real gray_subcfi = gray_subbbc*gray_subji + gray_subbbi*gray_subji + gray_subbff*gray_subcah + gray_subbgj*gray_subcfd + gray_subbhc*gray_subcfd + gray_subbhf*gray_subcdg - gray_subbhh*gray_subcfd;
real gray_subcfj = gray_subbag + gray_subbce + gray_subbid*gray_subcfd + gray_subcbi - gray_subcee + gray_subcef + gray_subceg;
real gray_subcga = gray_subbda + gray_subbdg + gray_subbjc*gray_subcfd + gray_subcca - gray_subccb + gray_subccd + gray_subcej;
real16 GammaUPt = {-gray_subgf - gray_subic - gray_subjj,
gray_subbcf,
gray_subbdh,
gray_subbec,
gray_subbcf,
gray_subbdc*gray_subbgc + gray_subbeh*gray_subbej - gray_subbff*gray_subib - gray_subbfj*gray_subji + gray_subbgd*gray_subbge + gray_subbgh*gray_subbgi,
gray_subbhi,
gray_subbie,
gray_subbdh,
gray_subbhi,
-gray_subbdc*gray_subbhf + gray_subbej*gray_subbig + gray_subbge*gray_subbij - gray_subbih*gray_subji + gray_subbii*gray_subib + gray_subbja*gray_subbjb,
gray_subbjd,
gray_subbec,
gray_subbie,
gray_subbjd,
-gray_subbdc*gray_subbid + gray_subbge*gray_subbjh + gray_subbgi*gray_subcac + gray_subbjb*gray_subcaa - gray_subbjc*gray_subib + gray_subbje*gray_subji};
real16 GammaUPx = {-gray_subcae*gray_subhj - gray_subcaf*gray_subge - gray_subcah*gray_subjf,
gray_subcbc,
gray_subcbe,
gray_subcbf,
gray_subcbc,
gray_subbeh*gray_subcbg + gray_subbgc*gray_subcbd + gray_subbgd*gray_subbjb + gray_subbgh*gray_subcbj - gray_subcbh - gray_subcbi,
gray_subcce,
gray_subccj,
gray_subcbe,
gray_subcce,
-gray_subbhf*gray_subcbd + gray_subbig*gray_subcbg - gray_subbih*gray_subcah + gray_subbii*gray_subcae + gray_subbij*gray_subbjb + gray_subbja*gray_subcda,
gray_subcdb,
gray_subcbf,
gray_subccj,
gray_subcdb,
-gray_subbid*gray_subcbd + gray_subbjb*gray_subbjh - gray_subbjc*gray_subcae + gray_subbje*gray_subcah + gray_subcaa*gray_subcda + gray_subcac*gray_subcbj};
real16 GammaUPy = {-gray_subcad*gray_subge - gray_subcde*gray_subhj - gray_subcdg*gray_subjf,
gray_subcdh,
gray_subcea,
gray_subceb,
gray_subcdh,
gray_subbeh*gray_subcec - gray_subbff*gray_subcde - gray_subbfj*gray_subcdg + gray_subbgc*gray_subcae + gray_subbgd*gray_subbgi + gray_subbgh*gray_subced,
gray_subceh,
gray_subcei,
gray_subcea,
gray_subceh,
gray_subbgi*gray_subbij + gray_subbig*gray_subcec + gray_subbii*gray_subcde + gray_subbja*gray_subcbj - gray_subccc - gray_subcej,
gray_subcfb,
gray_subceb,
gray_subcei,
gray_subcfb,
gray_subbgi*gray_subbjh - gray_subbid*gray_subcae - gray_subbjc*gray_subcde + gray_subbje*gray_subcdg + gray_subcaa*gray_subcbj + gray_subcac*gray_subced};
real16 GammaUPz = {-gray_subcag*gray_subge - gray_subcdg*gray_subhj - gray_subcfd*gray_subjf,
gray_subcfe,
gray_subcff,
gray_subcfg,
gray_subcfe,
gray_subbeh*gray_subcfh + gray_subbej*gray_subbgd - gray_subbff*gray_subcdg - gray_subbfj*gray_subcfd + gray_subbgc*gray_subcah + gray_subbgh*gray_subcec,
gray_subcfi,
gray_subcfj,
gray_subcff,
gray_subcfi,
gray_subbej*gray_subbij - gray_subbhf*gray_subcah + gray_subbig*gray_subcfh - gray_subbih*gray_subcfd + gray_subbii*gray_subcdg + gray_subbja*gray_subcbg,
gray_subcga,
gray_subcfg,
gray_subcfj,
gray_subcga,
gray_subbej*gray_subbjh + gray_subbje*gray_subcfd + gray_subcaa*gray_subcbg + gray_subcac*gray_subcec - gray_subccf - gray_subcfa};
real4 rhs = {-dot(u, matrix_vector_product(GammaUPt, u)),-dot(u, matrix_vector_product(GammaUPx, u)),-dot(u, matrix_vector_product(GammaUPy, u)),-dot(u, matrix_vector_product(GammaUPz, u))};
