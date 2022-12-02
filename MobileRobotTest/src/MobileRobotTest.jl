module MobileRobot

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

    mutable struct RobotModel
        radius::Float64
        wheel_offset::Float64
    end

    mutable struct RobotState
        x::Array{Float64}
        dx::Array{Float64}
        x_history::Matrix{Float64}
        dx_history::Matrix{Float64}
        time_history::Array{Float64}
    end

    mutable struct Robot
        robot_model::RobotModel
        robot_state::RobotState
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
        G = [cos(state.x[3]) 0; sin(state.x[3]) 0; 0 1]
        state.dx = G * u
        state.x = state.x + state.dx * param.sampling_time
        state.x_history = vcat(state.x_history, reshape(state.x, (1, 3)))
        state.dx_history = vcat(state.dx_history, reshape(state.dx, (1, 3)))
        push!(state.time_history, time)
    end

    function controller(param::TestCondition, state::RobotState, k)

            # nominal controller
            v_des = [1.0, 1.0]

            return v_des
    end

    function simulate(param::TestCondition, gains, name)

        # initialize
        x = param.x_0
        dx = zeros(3)
        x_history = reshape(x, (1, 3))
        dx_history = reshape(dx, (1, 3))
        time_history = [0.0]
        state = RobotState(x, dx, x_history, dx_history, time_history)

        # start
        count = 0
        while true

            u = controller(param, state, gains...)

            # update robot state
            time = param.sampling_time * count
            update_state(param, state, u, time)

            # check goal state
            if check_goal(param, state) == true || count > 1000
                break
            end

            count = count + 1
        end

        # plot results
        plot_pic(param, state, name)
        plot_gif(param, state, name)

    end

    function circle_shape(x, y, r)
        theta = LinRange(0, 2pi, 500)
        x .+ r * sin.(theta), y .+ r * cos.(theta)
    end

    function plot_pic(param::TestCondition, state::RobotState, name)
        plot(state.x_history[:, 1], state.x_history[:, 2], label="trajectory", xlabel="x [m]", ylabel="y [m]")
        plot!(circle_shape(param.x_o1[1], param.x_o1[2], param.d_obs), seriestype=[:shape,], c=:blue, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 1")
        plot!(circle_shape(param.x_o2[1], param.x_o2[2], param.d_obs), seriestype=[:shape,], c=:purple, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 2")
        png(name)

        p1 = plot(state.time_history[:, 1], state.x_history[:, 1], ylabel="x [m]")
        p2 = plot(state.time_history[:, 1], state.x_history[:, 2], xlabel="time [s]", ylabel="y [m]")
        plot(p1, p2, layout = (2, 1), legend = false)
        png(name*"_pos_time")

        p3 = plot(state.time_history[:, 1], state.dx_history[:, 1], ylabel="dx [m/s]")
        p4 = plot(state.time_history[:, 1], state.dx_history[:, 2], xlabel="time [s]", ylabel="dy [m/s]")
        plot(p3, p4, layout = (2, 1), legend = false)
        png(name*"_vel_time")
    end

    function plot_gif(param::TestCondition, state::RobotState, name)
        # initialize a plot
        plt = plot(
            1,
            aspectratio=1,
            xlim=(-5,  5),
            ylim=(-5,  5),
            label="trajectory",
            xlabel="x [m]",
            ylabel="y [m]"
        )


        # build an animated gif by pushing new points to the plot
        anim = @animate for i=1:length(state.time_history)
            plot(plt, circle_shape(state.x_history[i, 1], state.x_history[i, 2], 0.1), seriestype=[:shape,], c=:blue, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Robot")
            push!(plt, state.x_history[i, 1], state.x_history[i, 2])

            # add obstacles
            plot!(circle_shape(param.x_o1[1], param.x_o1[2], param.d_obs), seriestype=[:shape,], c=:blue, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 1")
            plot!(circle_shape(param.x_o2[1], param.x_o2[2], param.d_obs), seriestype=[:shape,], c=:purple, linecolor=:black, fillalpha=0.2, aspectratio=1, label="Obstacle 2")

        end every 10
        gif(anim, name*"_xy_2d.gif", fps = 30)
    end

end # module

using .MobileRobot

# test condition
x_0 = [0.0, 0.0, 0.0]
x_goal = [3.0, 5.0, 0.0]
x_o1 = [1.0, 2.0]
x_o2 = [2.5, 3.0]
d_obs = 0.5
sampling_time = 0.01

param = MobileRobot.TestCondition(x_0, x_goal, x_o1, x_o2, d_obs, sampling_time)

# simulation
# parameters
k = 1.0
gains = [k]
MobileRobot.simulate(param, gains, "cbf")


