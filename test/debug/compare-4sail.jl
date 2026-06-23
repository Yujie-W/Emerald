# set some constants here
ks = 0.67;
ddb = 0.42;
ddf = 0.41;
sdb = 0.42;
sdf = 0.41;
rho = 0.3;
tau = 0.25;
sigb = ddb * rho + ddf * tau;
sigf = ddf * rho + ddb * tau;
sb = sdb * rho + sdf * tau;
sf = sdf * rho + sdb * tau;
att = 1.0 - sigf;
m = sqrt(att^2 - sigb^2);
rinf = (att - m) / sigb


function Jfunc1(k, l, t)
    del_ = (k - l) * t;
    @show del_;
    return abs(del_) > 1e-3 ? (exp(-l * t) - exp(-k * t)) / (k - l) : 0.5 * t * (exp(-k * t) + exp(-l * t)) * (1 - (del_^2) / 12)
end;

function Jfunc2(k, l, t)
    return (1 - exp(-(k + l) * t)) / (k + l)
end;





function rho_tau_dd(l)
    e = exp(-m * l);
    denom = 1.0 - rinf^2 * e^2;
    re = rinf * e;

    tdd = (1.0 - rinf^2) * e / denom;
    rdd = rinf * (1.0 - e^2) / denom;

    J1ks = Jfunc1(ks, m, l);
    J2ks = Jfunc2(ks, m, l);
    Pss = (sf + sb * rinf) * J1ks;
    Qss = (sf * rinf + sb) * J2ks;

    tsd = (Pss - re * Qss) / denom;
    rsd = (Qss - re * Pss) / denom;

    tss = exp(-ks * l);

    return tdd, rdd, tsd, rsd, tss
end;

function rho_tau_double(l)
    # double adding algorithm with ndb = 10
    r_sd = 0;
    t_sd = 0;
    r_dd = 0;
    t_dd = 0;
    flai_10 = FT(2 ^ -20);
    r_dd = sigb * flai_10 * l;
    t_dd = sigf * flai_10 * l + 1 - flai_10 * l;
    a_dd = 1 - r_dd - t_dd;
    r_sd = sb * flai_10 * ks * l;
    t_sd = sf * flai_10 * ks * l;
    t_ss =  1 - flai_10 * ks * l;
    a_ss = 1 - r_sd - t_sd - t_ss;
    for idb in 1:20
        r_dd_2 = r_dd + t_dd * r_dd * t_dd / (1 - r_dd * r_dd);
        t_dd_2 = t_dd * t_dd / (1 - r_dd * r_dd);
        a_dd_2 = a_dd + t_dd * (1 + r_dd) * a_dd / (1 - r_dd * r_dd);
        @info "DD" r_dd_2 + t_dd_2 + a_dd_2;
        r_sd_2 = r_sd + t_sd * r_dd * t_dd / (1 - r_dd * r_dd) + t_ss * r_sd * t_dd / (1 - r_dd * r_dd);
        t_sd_2 = t_sd * t_dd / (1 - r_dd * r_dd) + t_ss * t_sd + t_ss * r_sd * r_dd * t_dd / (1 - r_dd * r_dd);
        t_ss_2 = t_ss * t_ss;
        a_ss_2 = a_ss + t_sd * (1 + r_dd) * a_dd / (1 - r_dd * r_dd) + t_ss * a_ss + t_ss * r_sd * (1 + r_dd) * a_dd / (1 - r_dd * r_dd);
        @info "SS" r_sd_2 + t_sd_2 + t_ss_2 + a_ss_2;
        r_dd = r_dd_2;
        t_dd = t_dd_2;
        a_dd = a_dd_2;
        r_sd = r_sd_2;
        t_sd = t_sd_2;
        t_ss = t_ss_2;
        a_ss = a_ss_2;
    end;

    return t_dd, r_dd, t_sd, r_sd, t_ss
end;

function rho_tau_2(tdd, rdd, tsd, rsd, tss)
    rdd2 = rdd + tdd * rdd * tdd / (1 - rdd * rdd)
    tdd2 = tdd * tdd / (1 - rdd * rdd)
    tsd2 = tsd * tdd / (1 - rdd * rdd) + tss * tsd + tss * rsd * rdd * tdd / (1 - rdd * rdd);
    rsd2 = rsd + tsd * rdd * tdd / (1 - rdd * rdd) + tss * rsd * tdd / (1 - rdd * rdd);

    return tdd2, rdd2, tsd2, rsd2
end;

for l in 0.5:0.5:3
    t1,r1,ts1,rs1,ss1 = rho_tau_dd(l);
    t2,r2,ts2,rs2,ss2 = rho_tau_dd(2l);
    t3,r3,ts3,rs3 = rho_tau_2(t1,r1,ts1,rs1,ss1);
    T1,R1,TS1,RS1,SS1 = rho_tau_double(l);
    T2,R2,TS2,RS2,SS2 = rho_tau_double(2l);
    T3,R3,TS3,RS3 = rho_tau_2(T1,R1,TS1,RS1,SS1);

    @info "debugging 1" l (t1,T1,t1/T1) (r1,R1,r1/R1) (ts1,TS1,ts1/TS1) (rs1,RS1,rs1/RS1);
    @info "debugging 2" l (t2,T2,t2/T2) (r2,R2,r2/R2) (ts2,TS2,ts2/TS2) (rs2,RS2,rs2/RS2);
    @info "debugging 3" l (t2,t3,t2/t3) (r2,r3,r2/r3) (ts2,ts3,ts2/ts3) (rs2,rs3,rs2/rs3);
    @info "debugging 4" l (T2,T3,T2/T3) (R2,R3,R2/R3) (TS2,TS3,TS2/TS3) (RS2,RS3,RS2/RS3);
end;
