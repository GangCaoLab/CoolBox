from coolbox.core.track.base import Track


class CustomTrack(Track):
    def __init__(self):
        super().__init__(properties_dict={})

    def fetch_data(self, gr, **kwargs):
        return "Test"

    def plot(self, ax, gr, **kwargs):
        x = (gr.start+gr.end)/2
        ax.text(x, 0, self.fetch_data(gr), fontsize=self.properties['fontsize'])
        ax.set_xlim(gr.start, gr.end)
