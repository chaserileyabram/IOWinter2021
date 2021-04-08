using Params.jl

# Chase Abram
# For solving HA discrete-time models via endogenous grid point



function solve_EGP(p, grids, income)

    @unpack_Params p
    # unpack other structs also?

    # Guess cons
    con_old = r.*aa_grid + yyF_grid + yyP_grid + yyT_grid
    p.con = con_old

    egp_iter = 0
    egp_diff = Inf

    while egp_iter < egp_maxiter && egp_diff > egp_tol

        # MU
        up(p)

        upinv(p)

        emuc(p)

        muc(p)

        update_con_euler(p)

        update_a_euler(p)

        interp_a_tom(p)

        update_con(p)

        egp_diff = maximum(abs.(p.con - con_old))

        con_old = p.con
        egp_iter += 1
    end





