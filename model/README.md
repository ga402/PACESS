```python
import torch

# Load trained YOLOv5 weights (best.pt)
model = torch.hub.load("ultralytics/yolov5", "custom", path="best.pt", source="local")

model.eval()  # inference mode
```
