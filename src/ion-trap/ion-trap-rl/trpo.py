import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, optimizers, losses, models
import numpy as np
# from utils import flatgrad, nn_model, assign_vars, flatvars
import os
import glob
from datetime import datetime
import threading
import gym
import time
import copy

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np

# Makes gradient of function loss_fn wrt var_list and
# flattens it to have a 1-D vector
def flatgrad(loss_fn, var_list):
    with tf.GradientTape() as t:
        loss = loss_fn()
    grads = t.gradient(loss, var_list, unconnected_gradients=tf.UnconnectedGradients.ZERO)
    return tf.concat([tf.reshape(g, [-1]) for g in grads], axis=0)

def nn_model(input_shape, output_shape, convolutional=False):
    model = keras.Sequential()
    if convolutional:
        model.add(layers.Lambda(lambda x: tf.cast(tf.image.resize(tf.image.rgb_to_grayscale(x), size=(32,32)), dtype=tf.float64)/256., input_shape=input_shape))
        model.add(layers.Conv2D(10, (3, 3), activation='relu'))
        model.add(layers.MaxPooling2D((3, 3)))
        model.add(layers.Conv2D(5, (3, 3), activation='relu'))
        model.add(layers.MaxPooling2D((3, 3)))
        model.add(layers.Flatten())
    # else:
    model.add(layers.Dense(64, input_shape=input_shape, activation='relu'))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(output_shape))
    return model

def nn_model2(input_shape, output_shape, convolutional=False):
    model = keras.Sequential()
    if convolutional:
        model.add(layers.Lambda(lambda x: tf.cast(tf.image.resize(tf.image.rgb_to_grayscale(tf.image.crop_to_bounding_box(x, 33,0,160,160)), size=(32,32)), dtype=tf.float64)/256., input_shape=input_shape))
        model.add(layers.Conv2D(20, (3, 3), activation='relu'))
        model.add(layers.MaxPooling2D((3, 3)))
        model.add(layers.Conv2D(20, (3, 3), activation='relu'))
        model.add(layers.MaxPooling2D((3, 3)))
        model.add(layers.Flatten())
    # else:
    model.add(layers.Dense(128, input_shape=input_shape, activation='relu'))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(output_shape))
    return model


def assign_vars(model, theta):
    """
    Create the process of assigning updated vars
    """
    shapes = [v.shape.as_list() for v in model.trainable_variables]
    size_theta = np.sum([np.prod(shape) for shape in shapes])
    # self.assign_weights_op = tf.assign(self.flat_weights, self.flat_wieghts_ph)
    start = 0
    for i, shape in enumerate(shapes):
        size = np.prod(shape)
        param = tf.reshape(theta[start:start + size], shape)
        model.trainable_variables[i].assign(param)
        start += size
    assert start == size_theta, "messy shapes"

def flatvars(model):
    return tf.concat([tf.reshape(v, [-1]) for v in model.trainable_variables], axis=0)

class TRPO:
    def __init__(self, env, training, policy_model=nn_model((15,), 1), value_model=nn_model((15,), 1), value_lr=1e-1, gamma=0.99, delta=0.01,
                 cg_damping=0.001, cg_iters=10, residual_tol=1e-5, ent_coeff=0.0, epsilon=0.4,
                 backtrack_coeff=0.6, backtrack_iters=10, batch_size=4096, epsilon_decay=lambda x: x - 5e-3, reward_scaling=1., correlated_epsilon=False):

        # current_time = datetime.now().strftime('%b%d_%H-%M-%S')
        # self.name = f"mylogs/TRPO-{self.env_name}-{current_time}"

        self.env = env

        self.epsilon_decay = epsilon_decay
        self.gamma = gamma
        self.cg_iters = cg_iters
        self.cg_damping = cg_damping
        self.ent_coeff = ent_coeff
        self.residual_tol = residual_tol
        self.delta = delta
        self.epsilon = epsilon
        self.backtrack_coeff = backtrack_coeff
        self.backtrack_iters = backtrack_iters
        self.reward_scaling = reward_scaling
        self.correlated_epsilon = correlated_epsilon
        self.BATCH_SIZE = batch_size

        self.model = policy_model
        self.tmp_model = models.clone_model(self.model)
        self.value_model = value_model

        if self.value_model:
            self.value_optimizer = optimizers.Adam(lr=value_lr)
            self.value_model.compile(self.value_optimizer, "mse")
            # self.writer = tf.summary.create_file_writer(self.name)

        self.training = training

        self.prev_state = self.env.reset()
        self.entropy = 0
        self.action = None
        self.action_prob = None
        self.action_list = []
        self.action_prob_list = []
        self.state_list = []
        self.reward_list = []
        self.G_list = []
        self.episodic_reward = 0.0

    def train_step(self, episode, obs_all, Gs_all, actions_all, action_probs_all, total_reward, best_reward, entropy, t0):

        def surrogate_loss(theta=None):
            if theta is None:
                model = self.model
            else:
                model = self.tmp_model
                assign_vars(self.tmp_model, theta)
            logits = model(obs)
            action_prob = tf.nn.softmax(logits)
            action_prob = tf.reduce_sum(actions_one_hot * action_prob, axis=1)
            old_logits = self.model(obs)
            old_action_prob = tf.nn.softmax(old_logits)
            old_action_prob = tf.reduce_sum(actions_one_hot * old_action_prob, axis=1).numpy() + 1e-8
            prob_ratio = action_prob / old_action_prob # pi(a|s) / pi_old(a|s)
            loss = tf.reduce_mean(prob_ratio * advantage) + self.ent_coeff * entropy
            return loss

        def kl_fn(theta=None):
            if theta is None:
                model = self.model
            else:
                model = self.tmp_model
                assign_vars(self.tmp_model, theta)
            logits = model(obs)
            action_prob = tf.nn.softmax(logits).numpy() + 1e-8
            old_logits = self.model(obs)
            old_action_prob = tf.nn.softmax(old_logits)
            return tf.reduce_mean(tf.reduce_sum(old_action_prob * tf.math.log(old_action_prob / action_prob), axis=1))

        def hessian_vector_product(p):
            def hvp_fn():
                kl_grad_vector = flatgrad(kl_fn, self.model.trainable_variables)
                grad_vector_product = tf.reduce_sum(kl_grad_vector * p)
                return grad_vector_product

            fisher_vector_product = flatgrad(hvp_fn, self.model.trainable_variables).numpy()
            return fisher_vector_product + (self.cg_damping * p)

        def conjugate_grad(Ax, b):
            """
            Conjugate gradient algorithm
            (see https://en.wikipedia.org/wiki/Conjugate_gradient_method)
            """
            x = np.zeros_like(b)
            r = b.copy() # Note: should be 'b - Ax(x)', but for x=0, Ax(x)=0. Change if doing warm start.
            p = r.copy()
            old_p = p.copy()
            r_dot_old = np.dot(r,r)
            for _ in range(self.cg_iters):
                z = Ax(p)
                alpha = r_dot_old / (np.dot(p, z) + 1e-8)
                old_x = x
                x += alpha * p
                r -= alpha * z
                r_dot_new = np.dot(r,r)
                beta = r_dot_new / (r_dot_old + 1e-8)
                r_dot_old = r_dot_new
                if r_dot_old < self.residual_tol:
                    break
                old_p = p.copy()
                p = r + beta * p
                if np.isnan(x).any():
                    print("x is nan")
                    print("z", np.isnan(z))
                    print("old_x", np.isnan(old_x))
                    print("kl_fn", np.isnan(kl_fn()))
            return x

        def linesearch(x, fullstep):
            fval = surrogate_loss(x)
            for (_n_backtracks, stepfrac) in enumerate(self.backtrack_coeff**np.arange(self.backtrack_iters)):
                xnew = x + stepfrac * fullstep
                newfval = surrogate_loss(xnew)
                kl_div = kl_fn(xnew)
                if np.isnan(kl_div):
                    print("kl is nan")
                    print("xnew", np.isnan(xnew))
                    print("x", np.isnan(x))
                    print("stepfrac", np.isnan(stepfrac))
                    print("fullstep",  np.isnan(fullstep))
                if kl_div <= self.delta and newfval >= 0:
                    print("Linesearch worked at ", _n_backtracks)
                    return xnew
                if _n_backtracks == self.backtrack_iters - 1:
                    print("Linesearch failed.", kl_div, newfval)
            return x

        NBATCHES = len(obs_all) // self.BATCH_SIZE
        if len(obs_all) < self.BATCH_SIZE:
            NBATCHES += 1

        for batch_id in range(NBATCHES):

            obs = obs_all[batch_id*self.BATCH_SIZE: (batch_id + 1)*self.BATCH_SIZE]
            Gs = Gs_all[batch_id*self.BATCH_SIZE: (batch_id + 1)*self.BATCH_SIZE]
            actions = actions_all[batch_id*self.BATCH_SIZE: (batch_id + 1)*self.BATCH_SIZE]
            action_probs = action_probs_all[batch_id*self.BATCH_SIZE: (batch_id + 1)*self.BATCH_SIZE]

            Vs = self.value_model(obs).numpy().flatten()
            # advantage = Gs
            advantage = Gs - Vs
            advantage = (advantage - advantage.mean())/(advantage.std() + 1e-8)
            actions_one_hot = tf.one_hot(actions, self.env.action_space.shape[0], dtype="float64")
            policy_loss = surrogate_loss()
            policy_gradient = flatgrad(surrogate_loss, self.model.trainable_variables).numpy()

            step_direction = conjugate_grad(hessian_vector_product, policy_gradient)

            shs = .5 * step_direction.dot(hessian_vector_product(step_direction).T)

            lm = np.sqrt(shs / self.delta) + 1e-8
            fullstep = step_direction / lm
            if np.isnan(fullstep).any():

                print("fullstep is nan")
                print("lm", lm)
                print("step_direction", step_direction)
                print("policy_gradient", policy_gradient)

            oldtheta = flatvars(self.model).numpy()

            theta = linesearch(oldtheta, fullstep)

            if np.isnan(theta).any():
                print("NaN detected. Skipping update...")
            else:
                assign_vars(self.model, theta)

            kl = kl_fn(oldtheta)

            history = self.value_model.fit(obs, Gs, epochs=5, verbose=0)
            value_loss = history.history["loss"][-1]

            print(f"Ep {episode}.{batch_id}: Rw_mean {total_reward} - Rw_best {best_reward} - PL {policy_loss} - VL {value_loss} - KL {kl} - epsilon {self.epsilon} - time {time.time() - t0}")

        # if self.value_model:
        #     writer = self.writer
        #     with writer.as_default():
        #         tf.summary.scalar("reward", total_reward, step=episode)
        #         tf.summary.scalar("best_reward", best_reward, step=episode)
        #         tf.summary.scalar("value_loss", value_loss, step=episode)
        #         tf.summary.scalar("policy_loss", policy_loss, step=episode)

        self.epsilon = self.epsilon_decay(self.epsilon)

    def take_action(self):

        prev_state = self.prev_state[np.newaxis, :]
        logits = self.model(prev_state)
        action_prob = tf.nn.softmax(logits).numpy().ravel()
        action = np.random.choice(range(action_prob.shape[0]), p=action_prob)

        # epsilon greedy
        if np.random.uniform(0, 1) < self.epsilon:
            if self.correlated_epsilon and np.random.uniform(0, 1) < 0.8 and self.action is not None:
                action = self.action
            else:
                action = np.random.randint(0, self.env.action_space.shape[0])

        self.action = action
        self.action_prob = action_prob

        print("printing action")
        print(self.action)
        print("printing action")

        return [self.action]

    def observe_and_learn(self, state, fidelity=None):

        assert self.value_model is not None

        _state, reward, _done, _info = self.env.step(state, fidelity)
        self.episodic_reward += reward

        if self.training:

            self.reward_list.append(reward / self.reward_scaling)
            self.state_list.append(self.prev_state)
            self.action_list.append(self.action)
            self.action_prob_list.append(self.action_prob)

            self.entropy += - tf.reduce_sum(self.action_prob * tf.math.log(self.action_prob))

            self.prev_state = state

            if fidelity is not None:

                G = 0
                for r in self.reward_list[::-1]:
                    G = r + self.gamma * G
                    self.G_list.insert(0, G)

                mean_total_reward = sum(self.reward_list)
                mean_entropy = np.mean(self.entropy / len(self.action_list))

                self.train_step(0, self.state_list, self.G_list, self.action_list, self.action_prob_list, mean_total_reward, mean_total_reward, mean_entropy, 0)

    def reset(self):

        self.prev_state = self.env.reset()
        self.entropy = 0
        self.action = None
        self.action_prob = None
        self.action_list = []
        self.action_prob_list = []
        self.state_list = []
        self.reward_list = []
        self.G_list = []
        self.episodic_reward = 0.0

    def save_weights(self, save_dir, termination_condition_status=False):  # todo: Needs correction

        #### Experimental ####

        # if episode % 10 == 0 and episode != 0 and self.value_model:
        # self.model.save_weights(f"{self.name}/{episode}.ckpt")

        ######################

        if (termination_condition_status):
            self.prev_actor_model.save_weights(save_dir + "actor.h5")
        else:
            self.actor_model.save_weights(save_dir + "actor.h5")

        # todo: Needs to comply with the "termination_condition"
        # self.critic_model.save_weights(save_dir + "critic.h5")

        # self.target_actor.save_weights(save_dir + "target_actor.h5")
        # self.target_critic.save_weights(save_dir + "target_critic.h5")

    def load_weights(self, path):
        self.model.load_weights(path)


    # def close(self):
    #
    #     for env in self.envs:
    #         env.close()

    # def render_episode(self, n=1):
    #
    #     for i in range(n):
    #         ob = self.envs[0].reset()
    #         done = False
    #         action = None
    #         while not done:
    #             self.envs[0].render()
    #             action, _ = self(ob, action)
    #             ob, r, done, info = self.envs[0].step(action)

    # def sample(self):
    #
    #     obs_all, actions_all, rs_all, action_probs_all, Gs_all = [None], [None], [None], [None], [None]
    #
    #     mean_total_reward = [None]
    #     mean_entropy = [None]
    #
    #     ## if len(glob.glob("render")) > 0:
    #     ##     self.render = True
    #     ## else:
    #     ##     self.render = False
    #
    #     ## if self.render:
    #     ##     self.render_episode()
    #
    #     #### #### #### ####
    #     entropy = 0
    #     obs, actions, rs, action_probs, Gs = [], [], [], [], []
    #     ob = self.env.reset()
    #     done = False
    #
    #     last_action = None
    #     while not done:
    #         action, action_prob = self(ob, last_action)
    #         new_ob, r, done, info = self.env.step(action)
    #         last_action = action
    #         rs.append(r/self.reward_scaling)
    #         obs.append(ob)
    #         actions.append(action)
    #         action_probs.append(action_prob)
    #         entropy += - tf.reduce_sum(action_prob*tf.math.log(action_prob))
    #         ob = new_ob
    #
    #     G = 0
    #     for r in rs[::-1]:
    #         G = r + self.gamma*G
    #         Gs.insert(0, G)
    #
    #     mean_total_reward[0] = sum(rs)
    #     entropy = entropy / len(actions)
    #     mean_entropy[0] = entropy
    #     obs_all[0] = obs
    #     actions_all[0] = actions
    #     rs_all[0] = rs
    #     action_probs_all[0] = action_probs
    #     Gs_all[0] = Gs
    #     #### #### #### ####
    #
    #     mean_entropy = np.mean(mean_entropy)
    #     best_reward = np.max(mean_total_reward)
    #     mean_total_reward = np.mean(mean_total_reward)
    #     Gs_all = np.concatenate(Gs_all)
    #     obs_all = np.concatenate(obs_all)
    #     rs_all = np.concatenate(rs_all)
    #     actions_all = np.concatenate(actions_all)
    #     action_probs_all = np.concatenate(action_probs_all)
    #
    #     return obs_all, Gs_all, mean_total_reward, best_reward, actions_all, action_probs_all, mean_entropy

    # def __call__(self, ob, last_action=None):
    #
    #     ob = ob[np.newaxis, :]
    #     logits = self.model(ob)
    #     action_prob = tf.nn.softmax(logits).numpy().ravel()
    #     action = np.random.choice(range(action_prob.shape[0]), p=action_prob)
    #     # epsilon greedy
    #     if np.random.uniform(0,1) < self.epsilon:
    #         if self.correlated_epsilon and np.random.uniform(0,1) < 0.8 and last_action is not None:
    #             action = last_action
    #         else:
    #             action = np.random.randint(0,self.env.action_space.n)
    #     self.last_action = action
    #
    #     return action, action_prob
