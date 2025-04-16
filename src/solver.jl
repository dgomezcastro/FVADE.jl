using ForwardDiff
using LinearSolve

function Newton(G, x0; abs_tol=1e-4, max_iters=1e2)
    DG(x) = ForwardDiff.jacobian(G, x)
    x = copy(x0)
    δ = 2 * abs_tol * ones(size(x))
    iter = 1

    while maximum(abs.(δ)) > abs_tol
        if iter == max_iters
            error("Newton did not converge")
        end
        δ = solve(LinearProblem(DG(x), G(x))).u # Solve DG(x) delta = G(x)
        x = x - δ
        @debug "iter=$iter, error = $(maximum(abs.(δ)))"
        iter = iter + 1
    end
    return x
end