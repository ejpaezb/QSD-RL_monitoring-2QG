"""
Title: Deep Deterministic Policy Gradient (DDPG)
Author: The DDPG implementation is adapted from (https://keras.io/examples/rl/ddpg_pendulum/)
Description: Implementing DDPG algorithm on the "Ion Trap" Problem.
"""

"""
## Introduction
**Deep Deterministic Policy Gradient (DDPG)** is a model-free off-policy algorithm for
learning continous actions.
It combines ideas from DPG (Deterministic Policy Gradient) and DQN (Deep Q-Network).
It uses Experience Replay and slow-learning target networks from DQN, and it is based on
DPG,
which can operate over continuous action spaces.
This tutorial closely follow this paper -
[Continuous control with deep reinforcement learning](https://arxiv.org/pdf/1509.02971.pdf)
## Problem
We are trying to solve the **Ion Trap** control problem.
In this setting, we can take only two actions: swing left or swing right.
What make this problem challenging for Q-Learning Algorithms is that actions
are **continuous** instead of being **discrete**. That is, instead of using two
discrete actions like `-1` or `+1`, we have to select from infinite actions
ranging from `-2` to `+2`.
## Quick theory
Just like the Actor-Critic method, we have two networks:
1. Actor - It proposes an action given a state.
2. Critic - It predicts if the action is good (positive value) or bad (negative value)
given a state and an action.
DDPG uses two more techniques not present in the original DQN:
**First, it uses two Target networks.**
**Why?** Because it add stability to training. In short, we are learning from estimated
targets and Target networks are updated slowly, hence keeping our estimated targets
stable.
Conceptually, this is like saying, "I have an idea of how to play this well,
I'm going to try it out for a bit until I find something better",
as opposed to saying "I'm going to re-learn how to play this entire game after every
move".
See this [StackOverflow answer](https://stackoverflow.com/a/54238556/13475679).
**Second, it uses Experience Replay.**
We store list of tuples `(state, action, reward, next_state)`, and instead of
learning only from recent experience, we learn from sampling all of our experience
accumulated so far.
Now, let's see how is it implemented.

"""
import gym
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras import models
import numpy as np

"""
We use [OpenAIGym](http://gym.openai.com/docs) to create the environment.
We will use the `upper_bound` parameter to scale our actions later.
"""

# problem = "Pendulum-v0"
# env = gym.make(problem)

"""
To implement better exploration by the Actor network, we use noisy perturbations,
specifically
an **Ornstein-Uhlenbeck process** for generating noise, as described in the paper.
It samples noise from a correlated normal distribution.
"""


class OUActionNoise:
    def __init__(self, mean, std_deviation, theta=0.15, dt=1e-2, x_initial=None):
        self.theta = theta
        self.mean = mean
        self.std_dev = std_deviation
        self.dt = dt
        self.x_initial = x_initial
        self.reset()

    def __call__(self):
        # Formula taken from https://www.wikipedia.org/wiki/Ornstein-Uhlenbeck_process.

        # x = (
        #     self.x_prev
        #     + self.theta * (self.mean - self.x_prev) * self.dt
        #     + self.std_dev * np.sqrt(self.dt) * np.random.normal(size=self.mean.shape)
        # )

        x = (
                0.0
                + self.theta * (self.mean - 0.0) * self.dt
                + self.std_dev * np.sqrt(self.dt) * np.random.normal(size=self.mean.shape)
        )

        # Store x into x_prev
        # Makes next noise dependent on current one
        self.x_prev = x
        return x

    def reset(self):
        if self.x_initial is not None:
            self.x_prev = self.x_initial
        else:
            self.x_prev = np.zeros_like(self.mean)


"""
The `Buffer` class implements Experience Replay.
---
![Algorithm](https://i.imgur.com/mS6iGyJ.jpg)
---
**Critic loss** - Mean Squared Error of `y - Q(s, a)`
where `y` is the expected return as seen by the Target network,
and `Q(s, a)` is action value predicted by the Critic network. `y` is a moving target
that the critic model tries to achieve; we make this target
stable by updating the Target model slowly.
**Actor loss** - This is computed using the mean of the value given by the Critic network
for the actions taken by the Actor network. We seek to maximize this quantity.
Hence we update the Actor network so that it produces actions that get
the maximum predicted value as seen by the Critic, for a given state.
"""

class DDPG:
    def __init__(
            self, env, training, critic_lr, actor_lr, gamma, tau, ou_noise_sd, buffer_capacity=100000, batch_size=64
    ):

        self.env = env
        self.counter = 0
        # self.pulse = np.array([-2.4111659, 1.53511, 1.21170287, -0.71186317, 1.20462246, 1.52825255, -2.42575168])
        # self.pulse = np.array([0.40457698, 1.05856782, 1.31974751, 1.05840049, 0.40443285])

        self.num_states = self.env.observation_space.shape[0]
        self.num_actions = self.env.action_space.shape[0]

        self.upper_bound = self.env.action_space.high[0]
        self.lower_bound = self.env.action_space.low[0]

        self.critic_lr = critic_lr
        self.actor_lr = actor_lr

        self.gamma = gamma  # Discount factor for future rewards
        self.tau = tau  # Used to update target networks

        self.actor_model = self.get_actor()
        self.critic_model = self.get_critic()

        self.target_actor = self.get_actor()
        self.target_critic = self.get_critic()

        self.critic_optimizer = tf.keras.optimizers.Adam(self.critic_lr)
        self.actor_optimizer = tf.keras.optimizers.Adam(self.actor_lr)

        # Making the weights equal initially
        self.target_actor.set_weights(self.actor_model.get_weights())
        self.target_critic.set_weights(self.critic_model.get_weights())

        self.training = training
        if self.training:
            self.ou_noise = OUActionNoise(mean=np.zeros(1), std_deviation=float(ou_noise_sd) * np.ones(1))
        else:
            self.ou_noise = None

        self.buffer_capacity = buffer_capacity  # Number of "experiences" to store at max
        self.batch_size = batch_size  # Num of tuples to train on.
        self.buffer_counter = 0  # Its tells us num of times record() was called.
        self.episodic_reward = 0

        # Instead of list of tuples as the exp.replay concept go
        # We use different np.arrays for each tuple element
        self.state_buffer = np.zeros((self.buffer_capacity, self.num_states))
        self.action_buffer = np.zeros((self.buffer_capacity, self.num_actions))
        self.reward_buffer = np.zeros((self.buffer_capacity, 1))
        self.next_state_buffer = np.zeros((self.buffer_capacity, self.num_states))

        self.prev_state = self.env.reset()
        self.action = 0.0
        self.episodic_reward = 0.0

        self.prev_actor_model = None
        self.temporary_buffer = []

    def policy(self, state):
        """
        `policy()` returns an action sampled from our Actor network plus some noise for
        exploration.
        """
        sampled_actions = tf.squeeze(self.actor_model(state))
        if self.ou_noise is not None:
            noise = np.asscalar(self.ou_noise())

            # Adding noise to action
            sampled_actions = sampled_actions.numpy() + noise
        else:
            sampled_actions = sampled_actions.numpy()

        # We make sure action is within bounds
        legal_action = np.clip(sampled_actions, self.lower_bound, self.upper_bound)

        return [np.squeeze(legal_action)]

    @tf.function
    def update_target(self):
        """
        This update target parameters slowly
        Based on rate `tau`, which is much less than one.
        """

        # update target actor
        for (a, b) in zip(self.target_actor.variables, self.actor_model.variables):
            a.assign(b * self.tau + a * (1 - self.tau))

        # update target critic
        for (a, b) in zip(self.target_critic.variables, self.critic_model.variables):
            a.assign(b * self.tau + a * (1 - self.tau))

    # Takes (s,a,r,s') obervation tuple as input   
    def record(self, obs_tuple):
        # Set index to zero if buffer_capacity is exceeded,
        # replacing old records
        index = self.buffer_counter % self.buffer_capacity

        self.state_buffer[index] = obs_tuple[0]
        self.action_buffer[index] = obs_tuple[1]
        self.reward_buffer[index] = obs_tuple[2]
        self.next_state_buffer[index] = obs_tuple[3]

        self.buffer_counter += 1

    # Eager execution is turned on by default in TensorFlow 2. Decorating with tf.function allows
    # TensorFlow to build a static graph out of the logic and computations in our function.
    # This provides a large speed up for blocks of code that contain many small TensorFlow operations such as this one.
    @tf.function
    def update(self, state_batch, action_batch, reward_batch, next_state_batch):
        # Training and updating Actor & Critic networks.
        # See Pseudo Code.
        with tf.GradientTape() as tape:
            target_actions = self.target_actor(next_state_batch, training=True)
            y = reward_batch + self.gamma * self.target_critic(
                [next_state_batch, target_actions], training=True
            )
            critic_value = self.critic_model([state_batch, action_batch], training=True)
            critic_loss = tf.math.reduce_mean(tf.math.square(y - critic_value))

        critic_grad = tape.gradient(critic_loss, self.critic_model.trainable_variables)
        self.critic_optimizer.apply_gradients(
            zip(critic_grad, self.critic_model.trainable_variables)
        )

        with tf.GradientTape() as tape:
            actions = self.actor_model(state_batch, training=True)
            critic_value = self.critic_model([state_batch, actions], training=True)
            # Used `-value` as we want to maximize the value given
            # by the critic for our actions
            actor_loss = -tf.math.reduce_mean(critic_value)

        actor_grad = tape.gradient(actor_loss, self.actor_model.trainable_variables)
        self.actor_optimizer.apply_gradients(
            zip(actor_grad, self.actor_model.trainable_variables)
        )

    # We compute the loss and update parameters
    def learn(self):
        # Get sampling range
        record_range = min(self.buffer_counter, self.buffer_capacity)
        # Randomly sample indices
        batch_indices = np.random.choice(record_range, self.batch_size)

        # Convert to tensors
        state_batch = tf.convert_to_tensor(self.state_buffer[batch_indices])
        action_batch = tf.convert_to_tensor(self.action_buffer[batch_indices])
        reward_batch = tf.convert_to_tensor(self.reward_buffer[batch_indices])
        reward_batch = tf.cast(reward_batch, dtype=tf.float32)
        next_state_batch = tf.convert_to_tensor(self.next_state_buffer[batch_indices])

        self.update(state_batch, action_batch, reward_batch, next_state_batch)

    def get_actor(self):
        """
        Here we define the Actor and Critic networks. These are basic Dense models
        with `ReLU` activation.
        Note: We need the initialization for last layer of the Actor to be between
        `-0.003` and `0.003` as this prevents us from getting `1` or `-1` output values in
        the initial stages, which would squash our gradients to zero,
        as we use the `tanh` activation.
        """

        # Initialize weights between -3e-3 and 3-e3
        last_init = tf.random_uniform_initializer(minval=-0.003, maxval=0.003)

        inputs = layers.Input(shape=(self.num_states,))
        out = layers.Dense(256, activation="relu")(inputs)
        out = layers.Dense(256, activation="relu")(out)
        outputs = layers.Dense(1, activation="tanh", kernel_initializer=last_init)(out)

        # Our upper bound is "max_rabi" for "Ion Trap".
        outputs = outputs * self.upper_bound
        model = tf.keras.Model(inputs, outputs)

        return model

    def get_critic(self):
        """
        Here we define the Actor and Critic networks. These are basic Dense models
        with `ReLU` activation.
        Note: We need the initialization for last layer of the Actor to be between
        `-0.003` and `0.003` as this prevents us from getting `1` or `-1` output values in
        the initial stages, which would squash our gradients to zero,
        as we use the `tanh` activation.
        """

        # State as input
        state_input = layers.Input(shape=(self.num_states))
        state_out = layers.Dense(16, activation="relu")(state_input)
        state_out = layers.Dense(32, activation="relu")(state_out)

        # Action as input
        action_input = layers.Input(shape=(self.num_actions))
        action_out = layers.Dense(32, activation="relu")(action_input)

        # Both are passed through seperate layer before concatenating
        concat = layers.Concatenate()([state_out, action_out])

        out = layers.Dense(256, activation="relu")(concat)
        out = layers.Dense(256, activation="relu")(out)
        outputs = layers.Dense(1)(out)

        # Outputs single value for give state-action
        model = tf.keras.Model([state_input, action_input], outputs)

        return model

    def take_action(self):

        tf_prev_state = tf.expand_dims(tf.convert_to_tensor(self.prev_state), 0)
        self.action = self.policy(tf_prev_state)
        # self.action = [self.pulse[self.counter]]
        # print(self.action)

        # print(self.num_actions)
        # if self.counter < 6:
        #     self.counter += 1
        # else:
        #     self.counter = 0

        return self.action

    def observe_and_learn(self, state, fidelity=None):
        
        _state, reward, _done, _info = self.env.step(state, fidelity)
        self.episodic_reward += reward

        if self.training:

            self.temporary_buffer.append((self.prev_state, self.action, reward, state))

            if fidelity is not None:
                self.prev_actor_model = models.clone_model(self.actor_model)
                self.prev_actor_model.set_weights(self.actor_model.get_weights())

                for record in self.temporary_buffer:
                    self.record(record)
                    self.learn()
                    self.update_target()

                self.temporary_buffer = []

        self.prev_state = state

    def reset(self):

        self.episodic_reward = 0.0
        self.prev_state = self.env.reset()

    def save_weights(self, save_dir, termination_condition_status=False):

        if (termination_condition_status):
            self.prev_actor_model.save_weights(save_dir + "actor.h5")
        else:
            self.actor_model.save_weights(save_dir + "actor.h5")

        # todo: Needs to comply with the "termination_condition"
        # self.critic_model.save_weights(save_dir + "critic.h5")

        # self.target_actor.save_weights(save_dir + "target_actor.h5")
        # self.target_critic.save_weights(save_dir + "target_critic.h5")

    def load_weights(self, save_dir):

        self.actor_model.load_weights(save_dir + "actor.h5")

        # todo: Needs to comply with the "termination_condition"
        # self.critic_model.load_weights(save_dir + "critic.h5")

        # self.target_actor.load_weights(save_dir + "target_actor.h5")
        # self.target_critic.load_weights(save_dir + "target_critic.h5")