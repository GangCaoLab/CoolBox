"""

`+` operation rules:

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

"""

from coolbox.core.track.base import Track
from coolbox.core.frame import Frame
from coolbox.core.feature import Feature
from coolbox.core.coverage.base import Coverage
from coolbox.core.browser import WidgetsPanel, Browser


def test_track_add_track():
    assert isinstance(Track({}) + Track({}), Frame)


def test_track_add_frame():
    assert isinstance(Frame() + Track({}), Frame)
    assert isinstance(Track({}) + Frame(), Frame)


def test_track_feature():
    assert isinstance(Track({}) + Feature(test=""), Track)
    assert isinstance(Feature(test="") + Track({}), Track)


def test_frame_feature():
    frame = Track({}) + Track({})
    assert isinstance(frame + Feature(test=""), Frame)
    assert isinstance(Feature(test="") + frame, Frame)


def test_track_coverage():
    assert isinstance(Track({}) + Coverage({}), Track)
    assert isinstance(Coverage({}) + Track({}), Track)


def test_frame_coverage():
    frame = Track({}) + Track({})
    assert isinstance(frame + Coverage({}), Frame)
    assert isinstance(Coverage({}) + frame, Frame)


def test_coverage_feature():
    assert isinstance(Coverage({}) + Feature(test=""), Coverage)
    assert isinstance(Feature(test="") + Coverage({}), Coverage)


def test_frame_widgetspanel():
    assert isinstance(Frame() + WidgetsPanel(), Browser)
    assert isinstance(WidgetsPanel() + Frame(), Browser)


def test_frame_compose():
    frame1 = Track({}) + Track({})
    frame2 = Track({}) + Track({})
    frame3 = Track({}) + Track({})
    frame = frame1 + frame2 + frame3
    assert len(frame.tracks) == 6
