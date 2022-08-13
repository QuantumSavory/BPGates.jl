function test_cnotperm()
    @testset "CNOTPerm" begin
        state = BellState(2)
        for _ in 1:10
            new_state = apply!(copy(state), rand(CNOTPerm,1,2))
            @test state == new_state
        end
    end
end

test_cnotperm()
