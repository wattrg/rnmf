-- geometry stuff
dim = 2
length = RealVec2(301.0, 301.0)
n_cells = UIntVec2(301, 301)

-- model
model = Model({
    H_far = RealVec2(1.0, 0.0),
    bubble_centre = RealVec2(301.0/2, 301.0/2),
    mu = RealVec2(3.2, 1.0),
    n_iter = 10,
    n_sub_iter = 550,
    relax = 0.2,
    bubble_radius = 301.0/20,
    tol = 0.1
})

-- actions
-- actions = Actions(
--     Action({name="solve", action = "solve", iterations=20})
-- )

