**Can I stop and restart `nomadic realtime`?**

Yes, you can safely stop and restart `nomadic realtime`. The program keeps track of what FASTQ files have been processed, and will continue processing from where it left off.

**Can I run more than one instance of `nomadic` at once?**

Yes, but if the process generates a webpage e.g. `nomadic dashboard` or `nomadic realtime`, only one will be live. You can run multiple instances of anything else e.g. `nomadic process` that doesn't have a dashboard.

