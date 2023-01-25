
a_state = Int(0b101011)

# a = BitArray(digits(a, base=2))
a = 9588514242
# a = 958851424237769245992348

function periodicity(a::Integer)

    a = digits(a, base=2)
    b = copy(a)
    c = copy(a)

    b = circshift(b,1)
    c = circshift(c,-1)
    n = 1

    while b != a && c != a
        b = circshift(b,1)
        c = circshift(c,-1)
        n += 1
    end
    return n
end

periodicity(a)



