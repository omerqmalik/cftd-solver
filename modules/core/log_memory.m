
function log_memory(from_time, to_time, message, varargin)
   global mem_log_file;
   persistent hostname;
   
   if isempty(hostname)
       [~,hostname] = system('hostname');
       hostname = regexprep(hostname,'\r\n|\n|\r','');
   end
   
   if isempty(mem_log_file)
       return
   end
   
   message = sprintf(message, varargin{:});
   %from_str = datestr(from,'yyyy-mm-dd,HH:MM:SS.fff');
   %to_str = datestr(to,'yyyy-mm-dd,HH:MM:SS.fff');
   log_file_id = fopen(mem_log_file,'at');
%    
%    
%    from_time.TimeZone = 'America/New_York';
%    to_time.TimeZone = 'America/New_York';
   
   total_mem = java.lang.Runtime.getRuntime.totalMemory;
   fprintf(log_file_id, '%s,%d,%10.0f,%10.0f,%s,%d\n', hostname, feature('getpid'), posixtime(from_time), posixtime(to_time), message, total_mem);
   fclose(log_file_id);
end

