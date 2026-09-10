**Can I stop and restart `nomadic realtime`?**

Yes, you can safely stop and restart `nomadic realtime`. The program keeps track of what FASTQ files have been processed, and will continue processing from where it left off.

**Can I run more than one instance of `nomadic` at once?**

Yes, just ensure you are not running nomadic multiple times using the same output folder, as those will interfere.
For `nomadic dashboard` or `nomadic realtime`, the dashboard will be opened using a free port.

**Can I use my own panel with `nomadic`?**

Yes, you can use your own panel with `nomadic` as long as it uses one of the currently supported reference genome. You will need to provide a bed file that describes the amplicon regions. Please refer to the [Advanced Usage](advanced.md) section for more details on how to configure `nomadic` with custom panels.
