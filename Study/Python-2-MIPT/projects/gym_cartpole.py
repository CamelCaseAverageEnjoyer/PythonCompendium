import gymnasium as gym
import numpy as np

env = gym.make("CartPole-v1", render_mode="rgb_array")
observation_space = env.observation_space.shape[0]
action_space = env.action_space.n   

terminal_state = False
state = env.reset()
# state = np.reshape(state, [1, observation_space])  

while True:
    env.render()
    # action = dqn_solver.act(state)
    action = np.random.randint(2)
    # print(env.step(action))
    state_next, reward, terminal, info, _ = env.step(action)
    state = np.reshape(state_next, [1, observation_space])
    # print(state)
      
    if terminal:
        break
