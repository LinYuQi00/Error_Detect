class Results:
    def __init__(self, point_id: str, err_val: float, err_abs: float):
        """
        初始化 Results 类
        :param point_id: 点的 ID
        :param err_val: 错误值
        :param err_abs: 错误的绝对值
        """
        self.point_id = point_id
        self.err_val = err_val
        self.err_abs = err_abs

    def __gt__(self, other):
        """
        重载 '>' 运算符，按错误绝对值比较大小
        :param other: 另一个 Results 对象
        :return: True 如果当前对象的 err_abs > other 的 err_abs
        """
        return self.err_abs > other.err_abs

    def __lt__(self, other):
        """
        重载 '<' 运算符，按 PointID 的数值大小比较
        :param other: 另一个 Results 对象
        :return: True 如果当前对象的 point_id 转为整数后 < other 的 point_id 转为整数后
        """
        return int(self.point_id) < int(other.point_id)

# 示例用法
if __name__ == "__main__":
    r1 = Results("101", 0.5, 1.2)
    r2 = Results("102", 0.3, 2.1)
    r3 = Results("99", 0.4, 0.9)

    # 比较 err_abs
    print(r1 > r2)  # False
    print(r3 < r1)  # True (按 PointID 的数值大小)
