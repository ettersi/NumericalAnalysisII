using Printf

# Utility functions

function pop_fast!(vec,idx)
    vec[idx],vec[end] = vec[end],vec[idx]
    return pop!(vec)
end

pop_rand!(vec) = pop_fast!(vec,rand(1:length(vec)))

function count_equals(board,i,j,di,dj)
    c = 1
    while i+c*di in 1:size(board,1) &&
          j+c*dj in 1:size(board,2) &&
          board[i+c*di,j+c*dj] == board[i,j]
        c += 1
    end
    return c - 1
end


"""
    winning_move(k,board,i,j)

Check whether the player who occupies `board[i,j]` has won the game by occupying
that square.
"""
winning_move(k,board,i,j) =
    count_equals(board,i,j, 0,-1) + count_equals(board,i,j, 0, 1) + 1 >= k ||
    count_equals(board,i,j,-1, 0) + count_equals(board,i,j, 1, 0) + 1 >= k ||
    count_equals(board,i,j,-1,-1) + count_equals(board,i,j, 1, 1) + 1 >= k ||
    count_equals(board,i,j,-1, 1) + count_equals(board,i,j, 1,-1) + 1 >= k


"""
    monte_carlo(m,n,k, n_samples) -> prob

Play `n_samples` rounds of the `m,n,k`-game to estimate the probabilities that
player 1 wins (`prob[1]`), player 2 wins (`prob[2]`), or the game ends in a draw
(`prob[3]`).
"""
function monte_carlo(m,n,k, n_samples)
    scores = zeros(3)
    for i = 1:n_samples
        winner = play_game(m,n,k)
        scores[winner] += 1
    end
    return scores./n_samples
end

"""
    play_game(m,n,k) -> winner

Play the `m,n,k`-game assuming random moves on behalf of both players,
"""
function play_game(m,n,k)
    board = zeros(Int,m,n)
    moves = vec([(i,j) for i = 1:m, j = 1:n])
    player = 1

    while !isempty(moves)
        i,j = pop_rand!(moves)
        board[i,j] = player
        if winning_move(k,board,i,j)
           return player
       end
       player = mod1(player+1,2)
    end

    return 3
end


"""
    recursion(m,n,k) -> prob

Compute the probabilities that the `m,n,k`-game ends with a win for player 1
(`prob[1]`), a win for player 2 (`prob[2]`), or a draw (`prob[3]`) by
enumerating all possible games.
"""
function recursion(m,n,k)
    board = zeros(Int,m,n)
    moves = vec([(i,j) for i = 1:m, j = 1:n])
    player = 1

    return recurse(k,board,moves,player)
end

function recurse(k,board,moves,player)
    n_moves = length(moves)
    if n_moves == 0
        # Base case: if we run out of moves, then the game must have ended in a draw.
        return [0.0,0.0,1.0]
    end

    prob = zeros(3)

    # Go through all possible moves in this turn
    for idx = 1:n_moves
        # Remove the current move from the list of possible moves, and update the board
        i,j = pop_fast!(moves,idx)
        board[i,j] = player

        if winning_move(k,board,i,j)
            prob[player] += 1/n_moves
        else
            prob .+= recurse(k,board,moves,mod1(player+1,2))./n_moves
        end

        # Reset the board and add the current move again to the list of possible move
        board[i,j] = 0
        push!(moves,(i,j))
        moves[idx],moves[end] = moves[end],moves[idx]
    end
    return prob
end



function main()
    m,n,k = 3,3,3
    n_samples = 100_000

    println("Runtimes:")
    time_mc  = @elapsed prob_mc  = monte_carlo(m,n,k, n_samples)
    println("  Monte Carlo: ", @sprintf("%.3f", time_mc), " sec")
    time_rec = @elapsed prob_rec = recursion(m,n,k)
    println("    Recursion: ", @sprintf("%.3f", time_rec), " sec")

    s = fill("",3,3)
    s[:,1] = (p->@sprintf("%.4f",p)).(prob_mc)
    s[:,2] = (p->@sprintf("%.4f",p)).(prob_rec)
    s[:,3] = (p->@sprintf("%.4f",p)).(abs.(prob_mc .- prob_rec))

    println()
    println("                   |   MC   | exact  | error  ")
    println("-------------------+--------+--------+--------")
    println("Win for player 1:  | ", s[1,1], " | ", s[1,2], " | ", s[1,3])
    println("Win for player 2:  | ", s[2,1], " | ", s[2,2], " | ", s[2,3])
    println("            Draw:  | ", s[3,1], " | ", s[3,2], " | ", s[3,3])
    println()
end
