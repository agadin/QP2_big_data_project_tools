import sys
import os
import math
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QWidget, QSlider, QFileDialog
from PyQt5.QtGui import QPixmap, QImage, QPainter, QPen
from PyQt5.QtCore import Qt, QPoint
from PIL import Image


class ImageMeasureApp(QMainWindow):
    def __init__(self, folder_path):
        super().__init__()
        self.folder_path = folder_path
        self.tif_files = self.get_tif_files(folder_path)
        self.current_index = 0
        self.init_ui()

    def get_tif_files(self, folder_path):
        """Retrieve all .tif files in the specified folder."""
        return sorted([f for f in os.listdir(folder_path) if f.endswith(".tif")])

    def init_ui(self):
        if not self.tif_files:
            print(f"No .tif files found in {self.folder_path}")
            sys.exit()

        self.load_image(self.tif_files[self.current_index])

        self.label = QLabel(self)
        self.label.setPixmap(self.pixmap)
        self.label.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.label.setScaledContents(True)
        self.label.mousePressEvent = self.get_point

        self.points = []

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(len(self.tif_files) - 1)
        self.slider.setValue(self.current_index)
        self.slider.valueChanged.connect(self.update_image)

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.slider)

        container = QWidget()
        container.setLayout(self.layout)
        self.setCentralWidget(container)

        self.setWindowTitle("Measure Distance on Image with File Navigation")
        self.resize(800, 800)
        self.show()

    def load_image(self, filename):
        """Load the image and prepare the pixmap."""
        image_path = os.path.join(self.folder_path, filename)
        self.image = Image.open(image_path)
        self.image_array = np.array(self.image)


        height, width = self.image_array.shape
        image_qt = QImage(
            self.image_array.data,
            width,
            height,
            self.image_array.strides[0],
            QImage.Format_Grayscale8
        )
        self.pixmap = QPixmap.fromImage(image_qt)
        self.original_pixmap = self.pixmap.copy()

    def update_image(self, value):
        """Update the displayed image when the slider value changes."""
        self.current_index = value
        self.load_image(self.tif_files[self.current_index])
        self.label.setPixmap(self.pixmap)
        self.points.clear()

    def get_point(self, event):
        """Capture click points on the image and calculate distance."""
        x = event.pos().x()
        y = event.pos().y()
        self.points.append(QPoint(x, y))

        if len(self.points) == 2:
            self.pixmap = self.original_pixmap.copy()
            painter = QPainter(self.pixmap)
            pen = QPen(Qt.blue, 3)
            painter.setPen(pen)
            painter.drawLine(self.points[0], self.points[1])
            painter.end()

            self.label.setPixmap(self.pixmap)
            x1, y1 = self.points[0].x(), self.points[0].y()
            x2, y2 = self.points[1].x(), self.points[1].y()
            pixel_distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            mm_distance = pixel_distance * (10 / 17.52)

            print(f"Distance: {pixel_distance:.2f} pixels ({mm_distance:.2f} mm)")
            self.points.clear()


# Main execution
if __name__ == "__main__":
    app = QApplication(sys.argv)

    # Specify the folder path containing .tif files
    folder_path = "./sample_group/CT_1/"
    window = ImageMeasureApp(folder_path)
    sys.exit(app.exec_())
