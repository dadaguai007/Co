 function dataout=BPF(datain,samplerate,LowF,HighF)
N=numel(datain);
data=reshape(datain,1,numel(datain));
t=1:1:N;
freq=samplerate*t/N;
filterResponse=zeros(size(freq));
position=find(freq<=HighF&freq>=LowF);
filterResponse(min(position):max(position))=1;
filterResponse(N-max(position)-1:N-min(position)-1)=1;
dataout=ifft(fft(data).*filterResponse);
dataout=reshape(dataout,size(datain));
end