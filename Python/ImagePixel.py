
class ImagePixel:
    def __init__(self, ID=-1, x=0.0, y=0.0):
        """
        Initializes an ImagePixel object.

        :param ID: Pixel ID.
        :param x: The x coordinate of the pixel.
        :param y: The y coordinate of the pixel.
        """
        self.ID = ID
        self.x = x
        self.y = y


# Example usage
if __name__ == '__main__':
    # Create an ImagePixel instance with default values
    default_pixel = ImagePixel()
    print(f"Default Pixel - ID: {default_pixel.ID}, x: {default_pixel.x}, y: {default_pixel.y}")

    # Create an ImagePixel instance with specific values
    specific_pixel = ImagePixel(ID=1, x=50.0, y=60.0)
    print(f"Specific Pixel - ID: {specific_pixel.ID}, x: {specific_pixel.x}, y: {specific_pixel.y}")
