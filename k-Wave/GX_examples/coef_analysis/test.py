def coef_train_v02_ld(gamma, w0, c, cfg, if_plot=False, message_step=1000):
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
def coef_train_v01_ld(gamma, w0, c, cfg, if_plot=False, message_step=1000):
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
