module CbfComparison

    using LinearAlgebra
    using Plots
    using Logging

    struct TestCondition
        x_0::Array{Float64}
        x_goal::Array{Float64}
        x_o1::Array{Float64}
        x_o2::Array{Float64}
        d_obs::Float64
        sampling_time::Float64
    end

    mutable struct RobotState
        x::Array{Float64}
        dx::Array{Float64}
        x_history::Matrix{Float64}
        dx_history::Matrix{Float64}
        time_history::Array{Float64}
    end

    function min_dist(param::TestCondition, state::RobotState)

        dist_o1 = norm(state.x - param.x_o1)
        dist_o2 = norm(state.x - param.x_o2)

        if dist_o1 < dist_o2
            dist = dist_o1 - param.d_obs
            return dist, state.x-param.x_o1
        else
            dist = dist_o2 - param.d_obs
            return dist, state.x-param.x_o2
        end

    end

    function check_goal(param::TestCondition, state::RobotState)

        dist_goal = norm(state.x - param.x_goal)

        if dist_goal < 0.01
            return true
        else
            return false
        end
    end

    function update_state(param::TestCondition, state::RobotState, u, time)

        state.dx = u
        state.x = state.x + state.dx * param.sampling_time
        state.x_history = vcat(state.x_history, reshape(state.x, (1, 2)))
        state.dx_history = vcat(state.dx_history, reshape(state.dx, (1, 2)))
        push!(state.time_history, time)

    end

    function apf(param::TestCondition, state::RobotState)

            # parameters
            k_att = 1.0
            k_rep = 1.0
            rho0 = 0.5

            # attractive potential
            delta_u_att = k_att * (state.x - param.x_goal)

            # repulsive potential
            rho, xor = min_dist(param, state)

            if rho > rho0
                delta_u_rep = zeros(2)
            else
                delta_u_rep = - k_rep/rho^2 * (1/rho - 1/rho0) * xor
            end

            # artificail potential
            f_apf = - delta_u_att - delta_u_rep

            return f_apf
    end

    function cbf()
        # place holder
    end

    function apf_cbf()
        # place holder
    end

    function simulate(param::TestCondition)

        # initialize
        x = param.x_0
        dx = zeros(2)
        x_history = reshape(x, (1, 2))
        dx_history = reshape(dx, (1, 2))
        time_history = [0.0]
        state = RobotState(x, dx, x_history, dx_history, time_history)

        # start
        count = 0
        while true

            # APF
            f_apf = apf(param, state)

            # update robot state
            time = param.sampling_time * count
            update_state(param, state, f_apf, time)

            # check goal state
            if check_goal(param, state) == true #|| count > 1000
                break  
            end

            count = count + 1
        end

        # plot results
        plot_pic(param, state)
        plot_gif(param, state)

    end

    function circle_shape(x, y, r)
        theta = LinRange(0, 2pi, 500)
        x .+ r * sin.(theta), y .+ r * cos.(theta)
    end

    function plot_pic(param::TestCondition, state::RobotState)
        plot(state.x_history[:, 1], state.x_history[:, 2], label="trajectory", xlabel="x [m]", ylabel="y [m]")
        plot!(circle_shape(param.x_o1[1], param.x_o1[2], param.d_obs), seriestype=[:shape,], c=:blue, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 1")
        plot!(circle_shape(param.x_o2[1], param.x_o2[2], param.d_obs), seriestype=[:shape,], c=:purple, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 2")
        png("xy_2d")

        p1 = plot(state.time_history[:, 1], state.x_history[:, 1], ylabel="x [m]")
        p2 = plot(state.time_history[:, 1], state.x_history[:, 2], xlabel="time [s]", ylabel="y [m]")
        plot(p1, p2, layout = (2, 1), legend = false)
        png("pos_time")

        p3 = plot(state.time_history[:, 1], state.dx_history[:, 1], ylabel="dx [m/s]")
        p4 = plot(state.time_history[:, 1], state.dx_history[:, 2], xlabel="time [s]", ylabel="dy [m/s]")
        plot(p3, p4, layout = (2, 1), legend = false)
        png("vel_time")
    end

    function plot_gif(param::TestCondition, state::RobotState)
        # initialize a plot
        plt = plot(
            1,
            label="trajectory", 
            xlabel="x [m]", 
            ylabel="y [m]"
        )

        # add obstacles
        plot!(circle_shape(param.x_o1[1], param.x_o1[2], param.d_obs), seriestype=[:shape,], c=:blue, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 1")
        plot!(circle_shape(param.x_o2[1], param.x_o2[2], param.d_obs), seriestype=[:shape,], c=:purple, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 2")

        # build an animated gif by pushing new points to the plot
        anim = @animate for i=1:length(state.time_history)
            push!(plt, state.x_history[i, 1], state.x_history[i, 2])
        end every 10
        gif(anim, "xy_2d.gif", fps = 30)
    end

end # module

using .CbfComparison

# test condition
x_0 = [0.0, 0.0]
x_goal = [3.0, 5.0]
x_o1 = [1.0, 2.0]
x_o2 = [2.5, 3.0]
d_obs = 0.5
sampling_time = 0.01

param = CbfComparison.TestCondition(x_0, x_goal, x_o1, x_o2, d_obs, sampling_time)

# start simulation
CbfComparison.simulate(param)