"""
CoolBox API

example:

``` Python
from coolbox.api import *

bsr = Browser()

# open files
bsr.open("GM12878-rna-1.bw")
bsr.open("GM12878-chip-h3k27ac-1.bw")
bsr.open("GM12878-Hind3-1.cool")

bsr.show()

```

Or, using another style(ggplot like) API.

``` Python
from coolbox.api import *

frame = XAxis() + \
        BigWig("GM12878-rna-1.bw) + Color("#66ccff") + \
        BigWig("GM12878-chip-h3k27ac-1.bw") + Color("#ff9c9c") + \
        Cool("GM12878-Hind3-1.cool") + ColorMap("#0000ff", "#ff0000")

bsr = Browser(frame)

bsr.show()

```


Element types:

* Track
* Coverage
* Frame
* Feature
    - FrameFeature
* WidgetsPanel
* Browser


The rule of element composition:
```
    Track + Track = Frame
    Track + Feature = Track
    Track + Coverage = Track
    Frame + Track = Frame
    Frame + Coverage = Frame
    Frame + Feature = Frame
    Frame + FrameFeature = Frame
    Frame + Frame = Frame
    Frame + WidgetsPanel = Browser
    Coverage + Feature = Coverage

    Frame * Feature = Frame
    Frame * Coverage = Frame
```
"""

import warnings

warnings.filterwarnings('ignore')

from coolbox.core.track import *
from coolbox.core.frame import *
from coolbox.core.feature import *
from coolbox.core.coverage import *
from coolbox.core.browser import *
from coolbox.core.frame.superframe import JointView
from coolbox.utilities.genome import GenomeRange
