import numpy as np
import math
from Image import *

class ErrorDetection:
    def __init__(self):
        self.m0 = 2.8e-3  # 毫米
        self.photographyScale = 12000
        self.img1 = Image()
        self.img2 = Image()
        self.cPixels = []  # 普通像素
        self.ePoints = set()  # 错误点
        self.GCPs = []  # 地面控制点

    def process(self, ioe_fileName, gcp_fileName, pixel_fileName, output_fileName):
        # Step 1: Read data
        self.read_test_file(ioe_fileName, gcp_fileName, pixel_fileName)

        # Step 2: Get common pixel points
        self.get_common_pixel()

        # Step 3: Perform relative orientation and error detection
        self.detect_error(output_fileName)

        return 0

    def read_test_file(self, ioe_file_name, gcp_file_name, pixel_file_name):
        """
        读取测试文件，包括内方位元素、地面控制点和像素数据。
        :param ioe_file_name: 内方位元素文件名
        :param gcp_file_name: 地面控制点文件名
        :param pixel_file_name: 像素数据文件名
        :return: 0 表示成功，-1 表示失败
        """
        # 读取内方位元素文件
        try:
            with open(ioe_file_name, 'r') as in_file:
                in_file.readline()  # 读取第一行（例如：1000 4）
                f, x0, y0 = map(float, in_file.readline().strip().split())
                inner_orientation = {"f": f, "x0": x0, "y0": y0}
        except FileNotFoundError:
            print("Open inner orientation elements file failed!")
            return -1

        # 读取地面控制点文件
        try:
            with open(gcp_file_name, 'r') as in_file:
                while True:
                    line = in_file.readline().strip()
                    if not line:
                        break
                    id_, y, x, z= map(float, line.split())
                    # id_, y, x, z, d = map(float, line.split())
                    if id_ == -99:
                        break
                    gcp = {"id": int(id_), "x": x, "y": y, "z": z}
                    # gcp = {"id": int(id_), "x": x, "y": y, "z": z, "d": int(d)}
                    # 控制点文件维度信息有则加上d
                    self.GCPs.append(gcp)
        except FileNotFoundError:
            print("Open ground control point file failed!")
            return -1

        # 读取像素文件
        try:
            with open(pixel_file_name, 'r') as in_file:
                # 读取第一幅图像像素
                in_file.readline()  # 跳过行
                id1 = int(in_file.readline().strip())
                id1 = abs(id1)
                in_file.readline()  # 跳过行
                while True:
                    line = in_file.readline().strip()
                    if not line:
                        break
                    ID, x, y = map(float, line.split())
                    if ID == -99:
                        break
                    self.img1.insert({"ID": int(ID), "x": x, "y": y})

                # 读取第二幅图像像素
                in_file.readline()  # 跳过行
                id2 = int(in_file.readline().strip())
                id2 = abs(id2)
                in_file.readline()  # 跳过行
                while True:
                    line = in_file.readline().strip()
                    if not line:
                        break
                    ID, x, y = map(float, line.split())
                    if ID == -99:
                        break
                    self.img2.add_pixel({"ID": int(ID), "x": x, "y": y})
        except FileNotFoundError:
            print("Open pixel file failed!")
            return -1

        # 更新图像信息
        self.img1.set_ioe(inner_orientation)
        self.img2.set_ioe(inner_orientation)
        self.img1.set_id(id1)
        self.img2.set_id(id2)

        return 0

    def get_common_pixel(self):
        """获取图像1和图像2中的公共像素"""
        for pixel in self.img1.pixels():
            if pixel in self.img2.pixels():
                self.cPixels.append(pixel)
        return 0

    def detect_error(self, output_filename):
        # Initialize variables
        num = len(self.cPixels)
        f = self.img1.get_inner().f  # Assuming this method exists to get 'f'
        phi, omega, kappa, u, v = 0, 0, 0, 0, 0
        A = np.zeros((num, 5))
        DX = np.zeros((5, 1))
        P = np.eye(num)  # Identity matrix for P
        Q = np.zeros((num, 1))

        loop_num1 = 0
        while True:
            loop_num1 += 1
            loop_num2 = 0
            while True:
                loop_num2 += 1

                # Calculate Bx, By, Bz
                Bx = 1
                By = Bx * u
                Bz = Bx * v
                a1 = math.cos(phi) * math.cos(kappa) - math.sin(phi) * math.sin(omega) * math.sin(kappa)
                a2 = -math.cos(phi) * math.sin(kappa) - math.sin(phi) * math.sin(omega) * math.cos(kappa)
                a3 = -math.sin(phi) * math.cos(omega)
                b1 = math.cos(omega) * math.sin(kappa)
                b2 = math.cos(omega) * math.cos(kappa)
                b3 = -math.sin(omega)
                c1 = math.sin(phi) * math.cos(kappa) + math.cos(phi) * math.sin(omega) * math.sin(kappa)
                c2 = -math.sin(phi) * math.sin(kappa) + math.cos(phi) * math.sin(omega) * math.cos(kappa)
                c3 = math.cos(phi) * math.cos(omega)

                # Process each point
                for i in range(num):
                    X1, Y1 = self.img1[self.cPixels[i]].x, self.img1[self.cPixels[i]].y
                    Z1 = -f
                    x2, y2 = self.img2[self.cPixels[i]].x, self.img2[self.cPixels[i]].y
                    X2 = a1 * x2 + a2 * y2 - a3 * f
                    Y2 = b1 * x2 + b2 * y2 - b3 * f
                    Z2 = c1 * x2 + c2 * y2 - c3 * f

                    # Point projection coefficients
                    N1 = (Bx * Z2 - Bz * X2) / (X1 * Z2 - Z1 * X2)
                    N2 = (Bx * Z1 - Bz * X1) / (X1 * Z2 - Z1 * X2)

                    # Fill matrix A and Q
                    A[i, 0] = -X2 * Y2 * N2 / Z2
                    A[i, 1] = -(Z2 + Y2 ** 2 / Z2) * N2
                    A[i, 2] = X2 * N2
                    A[i, 3] = Bx
                    A[i, 4] = -Y2 * Bx / Z2

                    Q[i, 0] = N1 * Y1 - N2 * Y2 - By

                # Solve the system of equations
                DX = np.linalg.inv(A.T @ P @ A) @ (A.T @ P @ Q)
                phi += DX[0, 0]
                omega += DX[1, 0]
                kappa += DX[2, 0]
                u += DX[3, 0]
                v += DX[4, 0]

                # Convergence check
                if np.max(np.abs(DX)) < 3e-5 or loop_num2 > 30:
                    break

            # Calculate residuals
            V = A @ DX - Q
            Qvv = np.linalg.inv(P) - A @ np.linalg.inv(A.T @ P @ A) @ A.T
            sigma0 = math.sqrt((V.T @ P @ V)[0, 0] / (num - 5))
            d = 3.5 + 82 / (81 + (sigma0 / 0.0028) ** 4)

            # Recalculate matrix P based on the method
            P_pre = P.copy()
            for i in range(num):
                Ti = V[i, 0] ** 2 / (sigma0 ** 2 * Qvv[i, i] * P[i, i])

                # Update P based on Ti
                wi = abs(V[i, 0]) / (sigma0 * math.sqrt(Qvv[i, i] * P[i, i]))

                if loop_num1 <= 3:
                    K = 1
                else:
                    K = 3.29

                # Update P matrix based on specific method conditions
                # Example method 1
                if Ti <= 0.01:
                    P[i, i] = P[i, i] + K * wi
                elif Ti <= 0.1:
                    P[i, i] = P[i, i] + 0.5 * K * wi
                else:
                    P[i, i] = P[i, i] + 0.1 * K * wi

                # If change in P is small, break the loop
                if np.max(np.abs(P_pre - P)) < 0.001 or loop_num2 > 30:
                    break

        # Calculate errors and write to output file
        error_results = []
        Bx = 1
        By = Bx * u
        Bz = Bx * v
        for point in self.ePoints:
            X1 = self.img1[point].x
            Y1 = self.img1[point].y
            Z1 = -f
            x2 = self.img2[point].x
            y2 = self.img2[point].y
            X2 = a1 * x2 + a2 * y2 - a3 * f
            Y2 = b1 * x2 + b2 * y2 - b3 * f
            Z2 = c1 * x2 + c2 * y2 - c3 * f

            N1 = (Bx * Z2 - Bz * X2) / (X1 * Z2 - Z1 * X2)
            N2 = (Bx * Z1 - Bz * X1) / (X1 * Z2 - Z1 * X2)

            dy = (N1 * Y1 - By) / N2 - Y2
            error_results.append({
                'PointID': str(point),
                'errAbs': abs(dy / self.m0),
                'errVal': dy / self.m0
            })

        # Write results to the output file
        with open(output_filename, 'w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['PointID', 'errAbs', 'errVal'])
            writer.writeheader()
            for result in error_results:
                writer.writerow(result)

        return 0







