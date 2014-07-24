years_tested = [2001 2004 2005 2006 2009 2012 2013];
constructs_tested = [1 3 15 12 32 547 1237];

years_snr = [2001 2005 2006 2009 2011 2012 2013];
snr = [0.4 0.5 3.75 7.5 2.5 12 30];

hf=figure(1);
set(hf, 'Position', [100 100 1000 500]);
hal = axes('Position', [.05 .1 .4 .8], 'FontSize', 16);
line(years_tested, constructs_tested);
har = axes('Position', [.55 .1 .4 .8], 'FontSize', 16);
line(years_snr, snr);

