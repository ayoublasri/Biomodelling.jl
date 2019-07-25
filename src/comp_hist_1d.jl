function comp_hist_1d(t, X, N)

# Compute the N+1 histogram of the process X [it starts from 0 up to N].
#
#  ---- INPUT PARAMETERS ----
#   t: time vector, [NoJx1]
#   X: state of the process as a function of time, [NoJx1]
#   N: maximum bin for the histogram, [1x1]
#
# ---- OUTPUT PARAMETERS ----
#   p_hist : histogram of frequencies of species population, [N+1x1]
#
# (c) Giorgos Arampatzis & Yannis Pantazis, 2012, umass & uoc
#     mail: {gtarabat, yannis.pantazis}@gmail.com

    NoJ = length(t)

    p_hist = zeros(N + 1, 1)

    dt = diff(t)

    for i = 0:N
        idx = X[1:NoJ - 1] == i
        p_hist[i + 1] = sum(dt[idx]) / t[end]
    end

    return p_hist

end