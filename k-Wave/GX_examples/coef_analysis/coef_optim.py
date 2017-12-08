import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt


def coef_train_v01(gamma, w0, c, n_x=100, lambda_ri=1., cutoff_perc=0, lr=3e-5, n_epoch=3000, if_plot=False, message_step=1000):
    # w^2 = c1 * k^2 + c2 * k^4 + c3 * (i*w) * k^2 + c4 * (i*w) * k^4

    ##################################################################
    # Generate the graph

    # Initialize input and coefficients to optimize
    x = tf.placeholder(dtype=tf.float32)
    k = tf.placeholder(dtype=tf.complex64)
    A1 = tf.Variable(1.0, dtype=tf.float32)
    A2 = tf.Variable(0.01, dtype=tf.float32)
    A3 = tf.Variable(0.01, dtype=tf.float32)
    A4 = tf.Variable(0.01, dtype=tf.float32)

    # Give the analytical solution of the equation
    e_a = tf.complex(A2, A4 * (1 + x))
    e_b = tf.complex(A1, A3 * (1 + x))
    e_c = tf.complex(-(1 + x) ** 2, 0.)
    delta = e_b ** 2 - 4 * e_a * e_c
    z = (-e_b + tf.sqrt(delta)) / (2 * e_a)

    # Derive the predicted phase velocity and attenuation coefficient
    zeta = z ** (1/2)
    k_pred = zeta * w0 / c
    w = (1 + x) * w0
    kr_pred = tf.real(k_pred)
    ki_pred = tf.imag(k_pred)
    cp_pred = w / kr_pred
    alpha_pred = -ki_pred

    # Define the loss function
    kr = tf.real(k)
    ki = tf.imag(k)
    cp = w / kr
    alpha = -ki
    loss = tf.reduce_sum(tf.square(cp_pred / cp - 1) + lambda_ri * tf.square(alpha_pred / alpha - 1))
    # loss = tf.reduce_sum(tf.square(cp - cp_pred) + lambda_ri * tf.square(alpha - alpha_pred))
    # loss = tf.reduce_sum(tf.square(kr - kr_pred) + lambda_ri * tf.square(ki - ki_pred))
    # optimizer = tf.train.GradientDescentOptimizer(lr)
    optimizer = tf.train.AdamOptimizer(lr)
    train = optimizer.minimize(loss)

    # Training data
    x_train = np.linspace(-1, 1, num=n_x)
    cutoff_b = np.ceil(x_train.shape[0] * cutoff_perc / 100)
    x_train = x_train[int(cutoff_b):]
    # x_train = x_train[1:]   # Remove the singularity at x = -1 (f = 0 Hz)

    k_train = w0 / c * (1 + x_train) ** (1 - gamma) * np.exp(1j * (-np.pi * gamma / 2))

    ##################################################################
    # Training loop

    A1_hist, A2_hist, A3_hist, A4_hist, loss_hist = [], [], [], [], []
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        for i_epoch in range(n_epoch):
            if i_epoch % message_step == 0:
                print('# epoch complete: ', i_epoch)
            A1_cur, A2_cur, A3_cur, A4_cur, loss_cur = sess.run([A1, A2, A3, A4, loss], {x: x_train, k: k_train})
            A1_hist.append(A1_cur)
            A2_hist.append(A2_cur)
            A3_hist.append(A3_cur)
            A4_hist.append(A4_cur)
            loss_hist.append(loss_cur)
            sess.run(train, {x: x_train, k: k_train})
        A1_con, A2_con, A3_con, A4_con, loss_con, cp_con, alpha_con, kr_con, ki_con, kr_pred_con, ki_pred_con = \
            sess.run([A1, A2, A3, A4, loss, cp_pred, alpha_pred, kr, ki, kr_pred, ki_pred], {x: x_train, k: k_train})

    ##################################################################
    # Validation

    print("A1: ", A1_hist[-1])
    print("A2: ", A2_hist[-1])
    print("A3: ", A3_hist[-1])
    print("A4: ", A4_hist[-1])
    print("loss: ", loss_hist[-1])

    f_axis = w0 * (1 + x_train) / (2 * np.pi)
    cp_the = w0 * (1 + x_train) / np.real(k_train)
    alpha_the = - np.imag(k_train)

    if if_plot:
        fig = plt.figure()
        plt.subplot(3, 2, 1)
        plt.plot(f_axis, kr_con, 'k', label='Theoretical')
        plt.plot(f_axis, kr_pred_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.subplot(3, 2, 2)
        plt.plot(f_axis, ki_con, 'k', label='Theoretical')
        plt.plot(f_axis, ki_pred_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.subplot(3, 2, 3)
        plt.plot(f_axis, cp_the, 'k', label='Theoretical')
        plt.plot(f_axis, cp_con, 'r--', label='Predicted')
        plt.subplot(3, 2, 4)
        plt.plot(f_axis, alpha_the, 'k', label='Theoretical')
        plt.plot(f_axis, alpha_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')
        plt.subplot(3, 2, 5)
        plt.plot(f_axis, (cp_con - cp_the) / cp_the * 100, 'r', label='Discrepancy (%)')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')
        plt.subplot(3, 2, 6)
        plt.plot(f_axis, (alpha_con - alpha_the) / alpha_the * 100, 'r', label='Discrepancy (%)')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')

        fig = plt.figure()
        plt.subplot(5, 1, 1)
        plt.plot(A1_hist)
        plt.ylabel('A1')
        plt.subplot(5, 1, 2)
        plt.plot(A2_hist)
        plt.ylabel('A2')
        plt.subplot(5, 1, 3)
        plt.plot(A3_hist)
        plt.ylabel('A3')
        plt.subplot(5, 1, 4)
        plt.plot(A4_hist)
        plt.ylabel('A4')
        plt.subplot(5, 1, 5)
        plt.plot(loss_hist)
        plt.ylabel('loss')

        plt.show()

    return

def coef_train_v02(gamma, w0, c, n_x=100, lambda_ri=1., cutoff_perc=0, lr=3e-5, n_epoch=3000, if_plot=False, message_step=1000):
    # w^2 = c1 * k + c2 * k^2 + c3 * (i*w) * k + c4 * (i*w) * k^2

    ##################################################################
    # Generate the graph

    # Initialize input and coefficients to optimize
    x = tf.placeholder(dtype=tf.float32)
    k = tf.placeholder(dtype=tf.complex64)
    A1 = tf.Variable(0.01, dtype=tf.float32)
    A2 = tf.Variable(1.00, dtype=tf.float32)
    A3 = tf.Variable(0.01, dtype=tf.float32)
    # A3 = tf.constant(0.00, dtype=tf.float32)
    A4 = tf.Variable(0.01, dtype=tf.float32)

    # Give the analytical solution of the equation
    e_a = tf.complex(A2, A4 * (1 + x))
    e_b = tf.complex(A1, A3 * (1 + x))
    e_c = tf.complex(-(1 + x) ** 2, 0.)
    delta = e_b ** 2 - 4 * e_a * e_c
    zeta = (-e_b + tf.sqrt(delta)) / (2 * e_a)

    # Derive the predicted phase velocity and attenuation coefficient
    k_pred = zeta * w0 / c
    w = (1 + x) * w0
    kr_pred = tf.real(k_pred)
    ki_pred = tf.imag(k_pred)
    cp_pred = w / kr_pred
    alpha_pred = -ki_pred

    # Define the loss function
    kr = tf.real(k)
    ki = tf.imag(k)
    cp = w / kr
    alpha = -ki
    loss = tf.reduce_sum(tf.square(cp_pred / cp - 1) + lambda_ri * tf.square(alpha_pred / alpha - 1))
    # loss = tf.reduce_sum(tf.square(cp - cp_pred) + lambda_ri * tf.square(alpha - alpha_pred))
    # loss = tf.reduce_sum(tf.square(kr - kr_pred) + lambda_ri * tf.square(ki - ki_pred))
    # optimizer = tf.train.GradientDescentOptimizer(lr)
    optimizer = tf.train.AdamOptimizer(lr)
    train = optimizer.minimize(loss)

    # Training data
    x_train = np.linspace(-1, 1, num=n_x)
    cutoff_b = np.ceil(x_train.shape[0] * cutoff_perc / 100)
    x_train = x_train[int(cutoff_b):]
    # x_train = x_train[1:]   # Remove the singularity at x = -1 (f = 0 Hz)
    k_train = w0 / c * (1 + x_train) ** (1 - gamma) * np.exp(1j * (-np.pi * gamma / 2))

    ##################################################################
    # Training loop

    A1_hist, A2_hist, A3_hist, A4_hist, loss_hist = [], [], [], [], []
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        for i_epoch in range(n_epoch):
            if i_epoch % message_step == 0:
                print('# epoch complete: ', i_epoch)
            A1_cur, A2_cur, A3_cur, A4_cur, loss_cur = sess.run([A1, A2, A3, A4, loss], {x: x_train, k: k_train})
            A1_hist.append(A1_cur)
            A2_hist.append(A2_cur)
            A3_hist.append(A3_cur)
            A4_hist.append(A4_cur)
            loss_hist.append(loss_cur)
            sess.run(train, {x: x_train, k: k_train})
        A1_con, A2_con, A3_con, A4_con, loss_con, cp_con, alpha_con, kr_con, ki_con, kr_pred_con, ki_pred_con = \
            sess.run([A1, A2, A3, A4, loss, cp_pred, alpha_pred, kr, ki, kr_pred, ki_pred], {x: x_train, k: k_train})

    ##################################################################
    # Validation

    print("A1: ", A1_hist[-1])
    print("A2: ", A2_hist[-1])
    print("A3: ", A3_hist[-1])
    print("A4: ", A4_hist[-1])
    print("loss: ", loss_hist[-1])

    f_axis = w0 * (1 + x_train) / (2 * np.pi)
    cp_the = w0 * (1 + x_train) / np.real(k_train)
    alpha_the = - np.imag(k_train)

    if if_plot:
        fig = plt.figure()
        plt.subplot(3, 2, 1)
        plt.plot(f_axis, kr_con, 'k', label='Theoretical')
        plt.plot(f_axis, kr_pred_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.subplot(3, 2, 2)
        plt.plot(f_axis, ki_con, 'k', label='Theoretical')
        plt.plot(f_axis, ki_pred_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.subplot(3, 2, 3)
        plt.plot(f_axis, cp_the, 'k', label='Theoretical')
        plt.plot(f_axis, cp_con, 'r--', label='Predicted')
        plt.subplot(3, 2, 4)
        plt.plot(f_axis, alpha_the, 'k', label='Theoretical')
        plt.plot(f_axis, alpha_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')
        plt.subplot(3, 2, 5)
        plt.plot(f_axis, (cp_con - cp_the) / cp_the * 100, 'r', label='Discrepancy (%)')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')
        plt.subplot(3, 2, 6)
        plt.plot(f_axis, (alpha_con - alpha_the) / alpha_the * 100, 'r', label='Discrepancy (%)')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')

        fig = plt.figure()
        plt.subplot(5, 1, 1)
        plt.plot(A1_hist)
        plt.ylabel('A1')
        plt.subplot(5, 1, 2)
        plt.plot(A2_hist)
        plt.ylabel('A2')
        plt.subplot(5, 1, 3)
        plt.plot(A3_hist)
        plt.ylabel('A3')
        plt.subplot(5, 1, 4)
        plt.plot(A4_hist)
        plt.ylabel('A4')
        plt.subplot(5, 1, 5)
        plt.plot(loss_hist)
        plt.ylabel('loss')

        plt.show()

    return A1_con, A2_con, A3_con, A4_con

def coef_train_v03(gamma, w0, c, n_x=100, lambda_ri=1., cutoff_perc=0, lr=3e-5, n_epoch=3000, if_plot=False, message_step=1000):
    # w^2 = c1 * k + c2 * k^2 + c3 * k^3 + c4 * (i*w) * k + c5 * (i*w) * k^2 + c6 * (i*w) * k^3

    ##################################################################
    # Generate the graph

    # Initialize input and coefficients to optimize
    x = tf.placeholder(dtype=tf.float32)
    k = tf.placeholder(dtype=tf.complex64)
    A1 = tf.Variable(-0.01, dtype=tf.float32)
    A2 = tf.Variable(1.00, dtype=tf.float32)
    A3 = tf.Variable(0.01, dtype=tf.float32)
    A4 = tf.Variable(0.01, dtype=tf.float32)
    A5 = tf.Variable(0.00, dtype=tf.float32)
    A6 = tf.Variable(0.00, dtype=tf.float32)

    # Give the analytical solution of the equation
    e_a = tf.complex(A3, A6 * (1 + x))
    e_b = tf.complex(A2, A5 * (1 + x))
    e_c = tf.complex(A1, A4 * (1 + x))
    e_d = tf.complex(-(1 + x) ** 2, 0.)

    delta0 = e_b ** 2 - 3 * e_a * e_c
    delta1 = 2 * e_b ** 3 - 9 * e_a * e_b * e_c + 27 * e_a ** 2 * e_d
    cc = ((delta1 + (delta1 ** 2 - 4 * delta0 ** 3) ** (1/2)) / 2) ** (1/3)
    # xi = tf.complex(1., 0.)
    # xi = tf.complex(-1 / 2, 1 / 2 * tf.sqrt(3.))
    xi = tf.complex(-1 / 2, -1 / 2 * tf.sqrt(3.))
    zeta = - (e_b + (xi * cc) + delta0 / (xi * cc)) / (3 * e_a)

    # Derive the predicted phase velocity and attenuation coefficient
    k_pred = zeta * w0 / c
    w = (1 + x) * w0
    kr_pred = tf.real(k_pred)
    ki_pred = tf.imag(k_pred)
    cp_pred = w / kr_pred
    alpha_pred = -ki_pred

    # Define the loss function
    kr = tf.real(k)
    ki = tf.imag(k)
    cp = w / kr
    alpha = -ki
    loss = tf.reduce_sum(tf.square(cp_pred / cp - 1) + lambda_ri * tf.square(alpha_pred / alpha - 1))
    # loss = tf.reduce_sum(tf.square(cp - cp_pred) + lambda_ri * tf.square(alpha - alpha_pred))
    # loss = tf.reduce_sum(tf.square(kr - kr_pred) + lambda_ri * tf.square(ki - ki_pred))
    # optimizer = tf.train.GradientDescentOptimizer(lr)
    optimizer = tf.train.AdamOptimizer(lr)
    train = optimizer.minimize(loss)

    # Training data
    x_train = np.linspace(-1, 1, num=n_x)
    cutoff_b = np.ceil(x_train.shape[0] * cutoff_perc / 100)
    x_train = x_train[int(cutoff_b):]
    # x_train = x_train[1:]   # Remove the singularity at x = -1 (f = 0 Hz)
    k_train = w0 / c * (1 + x_train) ** (1 - gamma) * np.exp(1j * (-np.pi * gamma / 2))

    ##################################################################
    # Training loop

    A1_hist, A2_hist, A3_hist, A4_hist, A5_hist, A6_hist, loss_hist = [], [], [], [], [], [], []
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        for i_epoch in range(n_epoch):
            if i_epoch % message_step == 0:
                print('# epoch complete: ', i_epoch)
            A1_cur, A2_cur, A3_cur, A4_cur, A5_cur, A6_cur, loss_cur = \
                sess.run([A1, A2, A3, A4, A5, A6, loss], {x: x_train, k: k_train})
            A1_hist.append(A1_cur)
            A2_hist.append(A2_cur)
            A3_hist.append(A3_cur)
            A4_hist.append(A4_cur)
            A5_hist.append(A5_cur)
            A6_hist.append(A6_cur)
            loss_hist.append(loss_cur)
            sess.run(train, {x: x_train, k: k_train})
        A1_con, A2_con, A3_con, A4_con, A5_con, A6_con, loss_con, \
        cp_con, alpha_con, kr_con, ki_con, kr_pred_con, ki_pred_con = \
            sess.run([A1, A2, A3, A4, A5, A6, loss, cp_pred, alpha_pred, kr, ki, kr_pred, ki_pred],
                     {x: x_train, k: k_train})

    ##################################################################
    # Validation

    print("A1: ", A1_con)
    print("A2: ", A2_con)
    print("A3: ", A3_con)
    print("A4: ", A4_con)
    print("A5: ", A5_con)
    print("A6: ", A6_con)
    print("loss: ", loss_con)

    f_axis = w0 * (1 + x_train) / (2 * np.pi)
    cp_the = w0 * (1 + x_train) / np.real(k_train)
    alpha_the = - np.imag(k_train)

    if if_plot:
        fig = plt.figure()
        plt.subplot(3, 2, 1)
        plt.plot(f_axis, kr_con, 'k', label='Theoretical')
        plt.plot(f_axis, kr_pred_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.subplot(3, 2, 2)
        plt.plot(f_axis, ki_con, 'k', label='Theoretical')
        plt.plot(f_axis, ki_pred_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.subplot(3, 2, 3)
        plt.plot(f_axis, cp_the, 'k', label='Theoretical')
        plt.plot(f_axis, cp_con, 'r--', label='Predicted')
        plt.subplot(3, 2, 4)
        plt.plot(f_axis, alpha_the, 'k', label='Theoretical')
        plt.plot(f_axis, alpha_con, 'r--', label='Predicted')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')
        plt.subplot(3, 2, 5)
        plt.plot(f_axis, (cp_con - cp_the) / cp_the * 100, 'r', label='Discrepancy (%)')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')
        plt.subplot(3, 2, 6)
        plt.plot(f_axis, (alpha_con - alpha_the) / alpha_the * 100, 'r', label='Discrepancy (%)')
        plt.legend(loc='best')
        plt.xlabel('Frequency (Hz)')

        fig = plt.figure()
        plt.subplot(7, 1, 1)
        plt.plot(A1_hist)
        plt.ylabel('A1')
        plt.subplot(7, 1, 2)
        plt.plot(A2_hist)
        plt.ylabel('A2')
        plt.subplot(7, 1, 3)
        plt.plot(A3_hist)
        plt.ylabel('A3')
        plt.subplot(7, 1, 4)
        plt.plot(A4_hist)
        plt.ylabel('A4')
        plt.subplot(7, 1, 5)
        plt.plot(A5_hist)
        plt.ylabel('A5')
        plt.subplot(7, 1, 6)
        plt.plot(A6_hist)
        plt.ylabel('A6')
        plt.subplot(7, 1, 7)
        plt.plot(loss_hist)
        plt.ylabel('loss')

        plt.show()

    return A1_con, A2_con, A3_con, A4_con, A5_con, A6_con

def main():

    # Set up the characteristic parameters of exploration seismology
    f0_m = 200.     # Reference frequency of the medium
    c0_m = 2089.    # Reference phase velocity
    rho  = 2200.    # Density
    Q = 32.5         # Quality value

    # Set up the frequency of interest
    f1 = 10.
    f2 = 60.
    f0 = 30.        # Reference frequency for calculation

    # Calculate the relevant parameters (Kjartansson's model)
    gamma = 1 / np.pi * np.arctan(1 / Q)
    w0 = 2 * np.pi * f0
    c0 = c0_m * (f0 / f0_m) ** gamma
    c = c0 * np.cos(np.pi * gamma / 2)

    # Hyper-Parmaeters
    lr = 1e-7           # Learning rate
    n_epoch = 20000      # Number of epochs (iterations)
    n_x = 100           # Number of samples
    lambda_ri = 1.      # Balance between real and imaginary parts
    cutoff_perc = 16

    # coef_train_v01(gamma, w0, c, n_x, lambda_ri, cutoff_perc, lr, n_epoch, if_plot=True)
    # coef_train_v02(gamma, w0, c, n_x, lambda_ri, cutoff_perc, lr, n_epoch, if_plot=True)
    coef_train_v03(gamma, w0, c, n_x, lambda_ri, cutoff_perc, lr, n_epoch, if_plot=True)

    # Loop over gamma list
    gamma_loop = False
    if gamma_loop:
        n_gamma = 10
        Q1, Q2 = 10., 100.
        gamma1 = 1 / np.pi * np.arctan(1 / Q1)
        gamma2 = 1 / np.pi * np.arctan(1 / Q2)
        gamma_vec = np.linspace(gamma1, gamma2, num=n_gamma)
        A1_vec = np.zeros(gamma_vec.shape)
        A2_vec = np.zeros(gamma_vec.shape)
        A3_vec = np.zeros(gamma_vec.shape)
        A4_vec = np.zeros(gamma_vec.shape)
        for i_gamma in range(n_gamma):
            print("Working on %d / %d" % (i_gamma + 1, n_gamma))
            gamma = gamma_vec[i_gamma]
            w0 = 2 * np.pi * f0
            c0 = c0_m * (f0 / f0_m) ** gamma
            c = c0 * np.cos(np.pi * gamma / 2)
            A1, A2, A3, A4 = coef_train_v02(gamma, w0, c, n_x, lambda_ri, cutoff_perc, lr, n_epoch)
            A1_vec[i_gamma] = A1
            A2_vec[i_gamma] = A2
            A3_vec[i_gamma] = A3
            A4_vec[i_gamma] = A4

        fig = plt.figure()
        plt.subplot(4, 1, 1)
        plt.plot(gamma_vec, A1_vec)
        plt.subplot(4, 1, 2)
        plt.plot(gamma_vec, A2_vec)
        plt.subplot(4, 1, 3)
        plt.plot(gamma_vec, A3_vec)
        plt.subplot(4, 1, 4)
        plt.plot(gamma_vec, A4_vec)
        plt.show()


    return

main()