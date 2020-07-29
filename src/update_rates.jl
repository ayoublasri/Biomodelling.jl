function update_rates(c, data, k_rea)

# the k-th flag is false if the k-th reaction rate is already updated.
    BOOL = ones(data.M, 1)

# loop over all species in "k_rea" reaction
    for i = 1:4

        spe_i = data.cm_rea[k_rea, i]

    # if i-th species is not zero
        if spe_i !== 0

        # loop over all reactions involving "spe_i"
            for k = 1:data.cm_spe[spe_i, 1]

                reaction = data.cm_spe[spe_i , k + 1]

            # if not already computed then compute rate
                if !iszero(BOOL[reaction])
                    c[reaction] =  rates(data, reaction)
                    BOOL[reaction] = 0
                end

            end # End "loop over all reactions involving 'spe_i'"

        end # End if

    end # End "loop over all species in 'k_rea' reaction"

    return c

end
