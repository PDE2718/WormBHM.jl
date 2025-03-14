function time_stamp()::String
    return string(now())
end
function rand_u32_tag(rng=TaskLocalRNG(),)::String
    return string(rand(rng, UInt32), base=16, pad=8)
end