ind = ~isnan(nLd(:,14)); %& abs(nLd(:,3))>2;
Redge = nLd(ind,end-2);
meanRedge = mean(Redge)
indb = nLd(ind,end);
meanindb = mean(indb)
the = mean(nLd(ind,end-1))

nLd(nLd(:,1)==0,:) = [];
nLd05 = nLd;
nLd1 = nLd;

ind05 = ~isnan(nLd05(:,14));% & abs(nLd05(:,3))<2;
Redge05 = nLd05(ind05,end-2);
meanRedge05 = mean(Redge05);
indb05 = nLd05(ind05,end);
meanindb05 = mean(indb05);
the05 = mean(nLd05(ind05,end-1));
meanthe05=mean(the05);
LC05 = nLd05(ind05,end-7);
meanLC05=mean(LC05);
Fcr05 = nLd05(ind05,5);
meanFcr05 = mean(Fcr05);

ind1 = ~isnan(nLd1(:,14));% & abs(nLd1(:,3))<2;
Redge1 = nLd1(ind1,end-2);
meanRedge1 = mean(Redge1)
indb1 = nLd1(ind1,end);
meanindb1 = mean(indb1);
the1 = mean(nLd1(ind1,end-1))
meanthe1=mean(the1);
LC1 = nLd1(ind1,end-7);
meanLC1=mean(LC1)
Fcr1 = nLd1(ind1,5);
meanFcr1 = mean(Fcr1);

quant05 = the05;
quant1 = the1;
(mean(quant05)-mean(quant1))/mean(quant05)
