## Some update constants
@kwdef struct UpdateConsts
    Cw::f64 = 0.5
    Eoff::f64 = 2.0
    P_del2ins::f64 = 1.0 # ratio p_delete/p_insert
end
# const Ecutoff::f64 = 1e3

@kwdef struct CycleAccumProb
    move_worm::f64   = 0.25
    insert_kink::f64 = 0.50
    delete_kink::f64 = 0.75
end
import Base.show
function Base.show(io::IO, a::CycleAccumProb)
    print(io, """Worm Cycle Probabilities
    P_move   = $(a.move_worm - 0.0)
    P_insert = $(a.insert_kink - a.move_worm)
    P_delete = $(a.delete_kink - a.insert_kink)
    P_glue   = $(1.0 - a.delete_kink)
    """
    )
end
function CycleProb(P_move_worm, P_insert_kink, P_delete_kink, P_glue_worm)::CycleAccumProb
    P_tot = P_move_worm + P_insert_kink + P_delete_kink + P_glue_worm
    p = CycleAccumProb(
        P_move_worm / P_tot,
        (P_move_worm + P_insert_kink) / P_tot,
        (P_move_worm + P_insert_kink + P_delete_kink) / P_tot,
    )
    @assert 0. < p.move_worm < p.insert_kink < p.delete_kink < 1.
    return p
end