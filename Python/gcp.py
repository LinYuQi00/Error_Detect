class GCP:
    class GCPpixel:
        def __init__(self, img_id=-1, x=0, y=0):
            self.img_id = img_id
            self.x = x
            self.y = y

    def __init__(self, id=None, x=None, y=None, z=None, demension=3):
        self.id = id if id is not None else -1
        self.x = x if x is not None else 0.0
        self.y = y if y is not None else 0.0
        self.z = z if z is not None else 0.0
        self.demension = demension
        self.pixels = []

    def get_id(self):
        return self.id

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_demension(self):
        return self.demension

    def add_pixel(self, img_id, x, y):
        """ Adds a GCPpixel to the GCP object """
        pixel = self.GCPpixel(img_id, x, y)
        self.pixels.append(pixel)

    def get_pixels(self):
        """ Returns the list of pixels associated with the GCP """
        return self.pixels
