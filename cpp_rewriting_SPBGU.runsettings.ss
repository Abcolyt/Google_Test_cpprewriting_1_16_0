<?xml version="1.0" encoding="utf-8"?>
<RunSettings>
  <!-- Конфигурация запуска -->
  <RunConfiguration>
    <ResultsDirectory>$(SolutionDir)TestResults</ResultsDirectory>
    <MaxCpuCount>1</MaxCpuCount>
    <EnvironmentVariables>
      <PATH>$(SolutionDir)..\vcpkg\installed\x64-windows\bin;$(PATH)</PATH>
    </EnvironmentVariables>
  </RunConfiguration>
  
  <!-- Настройки Google Test Adapter -->
  <GoogleTestAdapterSettings>
    <DebugMode>false</DebugMode>
    <ParallelTestExecution>true</ParallelTestExecution>
    <TraitsDiscovery>Enabled</TraitsDiscovery>
    
    <!-- Регулярное выражение для поиска тестовых исполняемых файлов -->
    <TestDiscoveryRegex>.*(Google_Test_cpprewriting_1_16_0|test).*\.exe</TestDiscoveryRegex>
    
    <!-- Дополнительные параметры -->
    <WorkingDir>$(ExecutableDir)</WorkingDir>
    <PathExtension>$(SolutionDir)..\vcpkg\installed\x64-windows\bin</PathExtension>
    <AdditionalTestExecutionParam>--gtest_output=xml</AdditionalTestExecutionParam>
    
    <!-- Таймауты -->
    <TestDiscoveryTimeoutInSeconds>30</TestDiscoveryTimeoutInSeconds>
    <TestExecutionTimeoutInSeconds>360</TestExecutionTimeoutInSeconds>
    
    <!-- Настройки фильтрации -->
    <SkipOriginCheck>false</SkipOriginCheck>
    <ShuffleTests>false</ShuffleTests>
    <ShuffleTestsSeed>0</ShuffleTestsSeed>
  </GoogleTestAdapterSettings>
  
  <!-- Настройки для Test Adapter for Google Test -->
  <TestAdapterForGoogleTestSettings>
    <DebugMode>false</DebugMode>
    <ParallelTestExecution>true</ParallelTestExecution>
    <TraitsDiscovery>Enabled</TraitsDiscovery>
    <TestDiscoveryTimeoutInSeconds>30</TestDiscoveryTimeoutInSeconds>
    <TestDiscoveryRegex>.*Google_Test_cpprewriting_1_16_0.*\.exe</TestDiscoveryRegex>
  </TestAdapterForGoogleTestSettings>
  
  <!-- Логирование -->
  <LoggerSettings>
    <LogLevel>0</LogLevel>
    <LogFilePath>$(SolutionDir)TestResults\test_log.txt</LogFilePath>
  </LoggerSettings>
  
  <!-- Конфигурации проектов -->
  <ProjectConfigurations>
    <Google_Test_cpprewriting_1_16_0>
      <Configuration>Debug</Configuration>
      <Platform>x64</Platform>
    </Google_Test_cpprewriting_1_16_0>
    <cpp_rewriting_SPBGU>
      <Configuration>Debug</Configuration>
      <Platform>x64</Platform>
    </cpp_rewriting_SPBGU>
  </ProjectConfigurations>
</RunSettings>