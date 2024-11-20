from ImagePixel import *

class Inner:
    def __init__(self, f=0.0, x0=0.0, y0=0.0):
        self.f = f
        self.x0 = x0
        self.y0 = y0

class Exterior:
    def __init__(self, Xs=0.0, Ys=0.0, Zs=0.0, phi=0.0, omega=0.0, kappa=0.0):
        self.Xs = Xs
        self.Ys = Ys
        self.Zs = Zs
        self.phi = phi
        self.omega = omega
        self.kappa = kappa

class Image:
    def __init__(self, ID=-1, in_data=None, ext_data=None, pixels=None):
        """
        Initializes the Image object.
        :param ID: The image ID.
        :param in_data: An instance of the Inner class.
        :param ext_data: An instance of the Exterior class.
        :param pixels: A dictionary mapping pixel IDs to ImagePixel objects.
        """
        if in_data is None:
            in_data = Inner()  # Use default Inner if None
        if ext_data is None:
            ext_data = Exterior()  # Use default Exterior if None
        if pixels is None:
            pixels = ImagePixel()

        self.ID = ID
        self.in_data = in_data
        self.ext_data = ext_data
        self.pixels = pixels
        self.pixelNum = len(pixels)

    # Change values
    def set_id(self, _ID):
        self.ID = _ID

    def set_ioe(self, in_data):
        self.in_data = in_data

    def set_eoe(self, ext_data):
        self.ext_data = ext_data

    # Get values
    def get_pixel_num(self):
        return self.pixelNum

    def get_id(self):
        return self.ID

    def get_inner(self):
        return self.in_data

    def get_ext(self):
        return self.ext_data

    # Access pixel by ID
    def get_pixel(self, _id):
        return self.pixels.get(_id)

    # Insert a new pixel
    def insert(self, pixel):
        self.pixels[pixel.ID] = pixel
        self.pixelNum += 1
        return 0

    # Find pixel by ID
    def find(self, _id):
        if _id in self.pixels:
            return _id
        else:
            return 0

# Example usage
if __name__ == '__main__':
    # Create Inner and Exterior instances
    inner = Inner(f=1.0, x0=5.0, y0=10.0)
    exterior = Exterior(Xs=10.0, Ys=20.0, Zs=30.0, phi=0.1, omega=0.2, kappa=0.3)

    # Create an Image instance
    image = Image(ID=123, in_data=inner, ext_data=exterior)

    # Create ImagePixel instance and insert into Image
    pixel = ImagePixel(ID=1, x=50, y=60)
    image.insert(pixel)

    # Retrieve and print pixel
    retrieved_pixel = image.get_pixel(1)
    print(f"Pixel ID: {retrieved_pixel.ID}, x: {retrieved_pixel.x}, y: {retrieved_pixel.y}")

    # Find pixel by ID
    found_pixel_id = image.find(1)
    print(f"Found pixel ID: {found_pixel_id}")
