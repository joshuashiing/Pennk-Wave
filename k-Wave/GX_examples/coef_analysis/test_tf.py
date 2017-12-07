import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

def main():

    # Hyper-parameters
    gamma = 0.01
    lr = 3e-5
    n_x = 100
    max_iter = 3000

    # Graph
    A1 = tf.Variable(1.0, dtype=tf.float32)
    A2 = tf.Variable(0.1, dtype=tf.float32)
    A3 = tf.Variable(0.1, dtype=tf.float32)
    A4 = tf.Variable(0.1, dtype=tf.float32)

    x = tf.placeholder(dtype=tf.float32)

    e_a = tf.complex(A2, A4 * (1 + x))
    e_b = tf.complex(A1, A3 * (1 + x))
    e_c = tf.complex(-(1 + x) ** 2, 0.)
    delta = e_b ** 2 - 4 * e_a * e_c
    z = (-e_b + tf.sqrt(delta)) / (2 * e_a)

    z_t = tf.placeholder(dtype=tf.complex64)
    loss = tf.reduce_sum(tf.square(tf.abs(z - z_t)))
    optimizer = tf.train.GradientDescentOptimizer(lr)
    train = optimizer.minimize(loss)

    # Training data
    x_train = np.linspace(-1, 1, num=n_x)
    z_train = ((1 + x_train) ** (2 - 2 * gamma)) * np.exp(1j * (-np.pi * gamma))

    # Training Loop
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    A1_hist = []
    A2_hist = []
    A3_hist = []
    A4_hist = []
    loss_hist = []
    for i in range(max_iter):
        sess.run(train, {x: x_train, z_t: z_train})
        curr_A1, curr_A2, curr_A3, curr_A4, curr_loss = \
            sess.run([A1, A2, A3, A4, loss], {x: x_train, z_t: z_train})
        A1_hist.append(curr_A1)
        A2_hist.append(curr_A2)
        A3_hist.append(curr_A3)
        A4_hist.append(curr_A4)
        loss_hist.append(curr_loss)

    print("A1: ", A1_hist[-1])
    print("A2: ", A2_hist[-1])
    print("A3: ", A3_hist[-1])
    print("A4: ", A4_hist[-1])
    print("loss: ", loss_hist[-1])

    z_solve = sess.run(z, {x: x_train})
    zeta_solve = np.sqrt(z_solve)
    zeta_train = np.sqrt(z_train)
    fig = plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(x_train, np.real(zeta_train), 'k')
    plt.plot(x_train, np.real(zeta_solve), 'r--')
    plt.subplot(2, 1, 2)
    plt.plot(x_train, np.imag(zeta_train), 'k')
    plt.plot(x_train, np.imag(zeta_solve), 'r--')
    plt.show()

    # fig = plt.figure()
    # plt.subplot(5, 1, 1)
    # plt.plot(A1_hist)
    # plt.subplot(5, 1, 2)
    # plt.plot(A2_hist)
    # plt.subplot(5, 1, 3)
    # plt.plot(A3_hist)
    # plt.subplot(5, 1, 4)
    # plt.plot(A4_hist)
    # plt.subplot(5, 1, 5)
    # plt.plot(loss_hist)
    # plt.show()





    return

main()