function comp_hist_1d(t, X, N)

    NoJ = length(t)

    p_hist = zeros(N + 1, 1)

    dt = diff(t)

    for i = 0:N
        idx = X[1:NoJ - 1] == i
        p_hist[i + 1] = sum(dt[idx]) / t[end]
    end

    return p_hist

end
