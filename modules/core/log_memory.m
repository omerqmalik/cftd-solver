
function log_memory(severity, message, varargin)
   global mem_log_file;
   if isempty(mem_log_file)
       return
   end
   
   message = sprintf(message, varargin{:});
   time_str = datestr(now,'yyyy-mm-dd HH:MM:SS.fff');
   log_file_id = fopen(mem_log_file,'at');
   
   total_mem = java.lang.Runtime.getRuntime.totalMemory;
   fprintf(log_file_id, '%s %d %s\t %s\n', time_str, total_mem, severity, message);
   fclose(log_file_id);
end

